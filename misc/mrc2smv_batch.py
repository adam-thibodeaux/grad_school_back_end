import os
import subprocess
input_dir = "/Users/adam/Documents/grad_school/danelius_lab/kobe_11152024_s1x1/raw_data/mrc"
output_dir = "/Users/adam/Documents/grad_school/danelius_lab/kobe_11152024_s1x1/raw_data"
sample_name = "s1x10"
count = 0
for f in sorted(os.listdir(input_dir)):
    if not os.path.isdir(f"{output_dir}/smv"):
        os.mkdir(f"{output_dir}/smv")
    mrc_string = f"mrc2smv -d 854 -r 1 -w '0.0251' -M 512 -E 1 -B 1 -o '{output_dir}/smv/{count}.img' {input_dir}/{f}"
    print(mrc_string)
    proc = subprocess.run([mrc_string],shell=True, capture_output=True)
    count += 1
    
