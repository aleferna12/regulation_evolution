runsdir=$1
outdir=$2
prefix=$3
framerate=$4

for rundir in "$runsdir"/*
do
  rundir=$(basename "$rundir")
	ffmpeg -framerate "$framerate" -pattern_type glob -i "$runsdir/$rundir/movie_/$prefix*.png" -c:v libx264 -sws_flags lanczos -pix_fmt yuv420p "$outdir/$rundir.mp4" > /dev/null &
done