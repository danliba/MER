# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 21:55:20 2023

@author: danli
"""

import cv2
import time

# Video path
video_path = "IsisROV_dive148_TRANSECT.m4v"

# Open the video
cap = cv2.VideoCapture(video_path)

# Check if video opened successfully
if not cap.isOpened():
    print("Error opening video file!")
    exit(1)

# Get video frame rate
fps = cap.get(cv2.CAP_PROP_FPS)

# Initialize counters
frame_count = 0
previous_time = time.time()  # Track time for interval messages

# Loop through frames
while True:
    # Read next frame
    success, frame = cap.read()

    if not success:
        break

    # Check if it's time to extract a frame and print interval message
    time_elapsed = time.time() - previous_time
    if frame_count == 0 or time_elapsed >= 10:
        previous_time = time.time()
        print(f"{(frame_count // fps) * 10} seconds, Frame number {frame_count + 1}")

    # Extract frame if needed
    if frame_count % fps == 0:
        # Create filename with frame count and timestamp
        filename = f"frame_{frame_count:05d}_{time.time():.2f}.jpg"

        # Save the frame as a JPEG image
        cv2.imwrite(filename, frame)
        print(f"Frame saved: {filename}")

    frame_count += 1

# Release the video capture
cap.release()

print("Go MER! Frame extraction finished!")

