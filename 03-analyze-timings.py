import glob
import re


def get_parameter(line, title):
    """Extracts the value corresponding to a specific title from a log file line."""
    if title not in line:
        return None
    parts = line.strip().split(": ", 1)
    return parts[1] if len(parts) > 1 else None


result = {}     # Dictionary to store results per offspring and method
methods = set()  # Set to store the names of all methods encountered

# Scan and parse all timing log files matching the pattern '*-timings*.log'
for filename in glob.glob('*-timings*.log'):
    method = filename.split('-')[0]  # Extract method name from filename
    methods.add(method)

    with open(filename, 'r') as file:
        offspring = ''
        time = 0
        hc_stage = 0  # Used specifically to handle HapCUT2's two-stage process
        stage1_time = 0
        stage1_memory = 0

        for line in file:
            value = get_parameter(line, 'Command being timed')
            if value:
                # For HapCUT2's extractHAIRS step
                if 'extractHAIRS' in value:
                    match = re.search(r'in/(\d+)_', value)
                    if match:
                        offspring = match.group(1)
                        hc_stage = 1  # First stage of HapCUT2
                # For HapCUT2's main phasing step
                elif 'HAPCUT2' in value:
                    if hc_stage == 1:
                        hc_stage += 1
                # For TinkerHap
                elif 'tinkerhap' in value:
                    match = re.search(r'in/(\d+)_', value)
                    if match:
                        offspring = match.group(1)
                # Fallback for other tools
                else:
                    match = re.search(r'out/\w*-(\d+)', value)
                    if match:
                        offspring = match.group(1)
                continue

            # Parse wall clock time
            value = get_parameter(line, 'wall clock')
            if value:
                parts = list(map(float, value.split(':')))
                time = sum(val * (60 ** i) for i, val in enumerate(reversed(parts)))
                continue

            # Parse memory usage
            value = get_parameter(line, 'Maximum resident set size')
            if value:
                memory = int(value)
                if method == 'HapCUT2':
                    if hc_stage == 1:
                        stage1_time = time
                        stage1_memory = memory
                    elif hc_stage == 2:
                        hc_stage = 0
                        result.setdefault(offspring, {})[method] = {
                            'time': time + stage1_time,
                            'memory': max(stage1_memory, memory)
                        }
                else:
                    result.setdefault(offspring, {})[method] = {'time': time, 'memory': memory}
                continue

# Collect all offspring IDs
offsprings = sorted(result.keys())
methods = sorted(methods)

# Prepare table headers
headers = ['sample']
for method in methods:
    headers.append(f'{method}-time')
    headers.append(f'{method}-mem')

output_lines = ['\t'.join(headers)]  # Header line for output file
times = {method: [] for method in methods}
mems = {method: [] for method in methods}

# Populate data per offspring
for offspring in offsprings:
    line = [offspring]
    for method in methods:
        time_val = result.get(offspring, {}).get(method, {}).get('time', '')
        mem_val = result.get(offspring, {}).get(method, {}).get('memory', '')
        if time_val != '':
            times[method].append(time_val)
        if mem_val != '':
            mems[method].append(mem_val)
        # Format time as float with two decimals
        rounded_time = f"{time_val:.2f}".rstrip('0').rstrip('.') if time_val != '' else ''
        line.extend([rounded_time, str(mem_val)])
    output_lines.append('\t'.join(line))

# Write results to output file
with open('output2.txt', 'w') as f:
    f.write('\n'.join(output_lines) + '\n')

# Print summary statistics
print(f"Parsed {len(offsprings)} samples")

for method in methods:
    if times[method]:
        times[method].sort()
        mems[method].sort()
        median_index = len(times[method]) // 2
        median_time = times[method][median_index]
        median_mem = mems[method][median_index]
        print(f"  {method}: Count={len(times[method])}; Time={median_time}; Mem={int(median_mem / 1024):,}")

        # List of samples missing for this method (if needed)
        missing = [offspring for offspring in offsprings if method not in result.get(offspring, {})]
