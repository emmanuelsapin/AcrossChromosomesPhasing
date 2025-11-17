#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Read the original file
with open('ProgramPhasing.cpp', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find the start and end of predict_phasing_without_parents (the full implementation, not the wrapper)
start_line = None
end_line = None
for i, line in enumerate(lines):
    if line.strip().startswith('int predict_phasing_without_parents(') and 'ID,int chr1' in line:
        start_line = i
    if start_line is not None and line.strip().startswith('int predict_phasing_with_parents_providing_GT(') and 'ID,int chr1' in line:
        end_line = i
        break

if start_line is None or end_line is None:
    print(f"Error: Could not find function boundaries. start={start_line}, end={end_line}")
    exit(1)

print(f"Extracting predict_phasing_without_parents from line {start_line+1} to {end_line}")

# Extract function 1
func1_lines = lines[start_line:end_line]
func1_code = ''.join(func1_lines)

# Adapt the code: replace global variable references with parameter references
# IDjob -> *IDjob (it's passed as int*)
func1_code = func1_code.replace('IDjob', '*IDjob')
# placefirttoconsider -> *placefirttoconsider (it's passed as int*)
func1_code = func1_code.replace('placefirttoconsider', '*placefirttoconsider')

# Write to file
with open('temp_func1.txt', 'w', encoding='utf-8') as f:
    f.write(func1_code)

# Find the start and end of predict_phasing_with_parents_providing_GT (the full implementation)
start_line2 = None
end_line2 = None
for i, line in enumerate(lines):
    if line.strip().startswith('int predict_phasing_with_parents_providing_GT(') and 'ID,int chr1' in line:
        start_line2 = i
    if start_line2 is not None and i > start_line2 and line.strip().startswith('int ') and '(' in line and ')' in line and not 'predict_phasing' in line:
        # Check if it's a function definition (not just a variable declaration)
        if '{' in lines[i+1] if i+1 < len(lines) else False:
            end_line2 = i
            break

if end_line2 is None:
    end_line2 = len(lines) - 1

print(f"Extracting predict_phasing_with_parents_providing_GT from line {start_line2+1} to {end_line2+1}")

# Extract function 2
func2_lines = lines[start_line2:end_line2+1]
func2_code = ''.join(func2_lines)

# Adapt the code
func2_code = func2_code.replace('IDjob', '*IDjob')
func2_code = func2_code.replace('placefirttoconsider', '*placefirttoconsider')

# Write to file
with open('temp_func2.txt', 'w', encoding='utf-8') as f:
    f.write(func2_code)

print("Extraction complete!")


