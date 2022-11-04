module purge
module load beta-modules
module load r/r-4.1.1_bc-3.13
module load nano/2.4.2
module load gerun
module load grass/8.0dev

export TZ="Europe/London"
export TMPDIR="/home/tcrnbgh/Scratch/tmp"

echo "TMPDIR=/home/tcrnbgh/Scratch/tmp" > /home/tcrnbgh/Scratch/visual_prominence/.Renviron
cp /home/tcrnbgh/Scratch/visual_prominence/.Renviron /home/tcrnbgh/.Renviron
chmod +x /home/tcrnbgh/Scratch/visual_prominence/python/viz.py



