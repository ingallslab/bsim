#!/bin/bash
#This script converts a series of images into an mp4 video
#ffmpeg -r fps -f image2 -s resolution (1920x1080) -i pic%04d.png -vcodec libx264 -crf quality (15-25)  -pix_fmt pixel_format videoName.mp4

cd "results/"

#rename all images in the folder so that ffmpeg can read them
num=1
for file in *.png; do
	printf -v n "%05d" $num
	mv "$file" "$n.png"
	let num=$num+1
done

ffmpeg -framerate 10 -s 640x360 -i %05d.png -codec:v libx264 test.mp4
