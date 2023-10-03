#!/bin/sh

docker container stop corems-molformdb-1
docker container rm corems-molformdb-1
docker volume rm corems_db-volume
docker compose -f /Users/christiandewey/CoreMS/docker-compose.yml up -d

