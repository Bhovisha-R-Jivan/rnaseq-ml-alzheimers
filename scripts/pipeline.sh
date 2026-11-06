# Step 1 = checking the quality of the raw reads using FastQC
bash /home/bhovs/rnaseq-ml-alzheimers/scripts/2.quality_assessment.sh

# Step 2 = trimming the adapters and low-quality bases using Trim Galore
bash /home/bhovs/rnaseq-ml-alzheimers/scripts/3.adapter_trimming.sh