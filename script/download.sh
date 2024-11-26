#!/bin/bash
#SBATCH --job-name=download_job      # Job name
#SBATCH --output=download_job_%j.log # Output log file (%j will be replaced with job ID)
#SBATCH --error=download_job_%j.err  # Error log file (%j will be replaced with job ID)
#SBATCH --ntasks=1                   # Number of tasks (processes)
#SBATCH --time=00:30:00              # Time limit (hh:mm:ss)
#SBATCH --mem=1G                     # Memory per node (adjust as needed)
#SBATCH --partition=standard         # Partition to submit to (adjust to your cluster setup)

# Load wget
module load wget

# Text file containing URLs
LINK_FILE="/scratch_tmp/grp/msc_appbio/practice/links2download.txt"

# Check if the file exists
if [[ ! -f "$LINK_FILE" ]]; then
echo "Error: $LINK_FILE does not exist."
exit 1
fi

# Create a directory to store the downloaded files
DOWNLOAD_DIR="/scratch_tmp/grp/msc_appbio/practice/downloads"
mkdir -p "$DOWNLOAD_DIR"

# Read each line of URL from the file and download it
while IFS= read -r URL; do
if [[ -n "$URL" ]]; then
echo "Downloading: $URL"
wget -P "$DOWNLOAD_DIR" "$URL" || echo "Failed to download: $URL"
fi
done < "$LINK_FILE"

echo "Download complete. Files saved in '$DOWNLOAD_DIR'."

