#!/bin/bash

#build --no-cache --pull
docker build -t msbase .
# run
docker run -it \
   -p 11434:11434 \
   -w /app \
   --memory=6g \
   --pids-limit=500 \
   msbase
