# Render video of lens.f90 from Introduction to Conventional Transmission Electron Microscopy
- https://ctem.web.cmu.edu
- https://ctem.web.cmu.edu/frames.html

## Instructions
- edit main.py and edit the first_frame and last_frame parameter in the main function
- the output video will interpolate between these parameters
- run `sudo docker build -t chrisbesch/ctemsoft . && sudo docker run -v ./out:/python_src/out chrisbesch/ctemsoft`
- the output will be in `out/lens.mp4`
