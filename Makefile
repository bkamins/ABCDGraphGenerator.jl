.PHONY: run
run:
	docker-compose run --rm app ash

.PHONY: save
save:
	docker save julia_abcd | gzip -9 > julia_abcd.tar.gz

.PHONY: load
load:
	docker load < ./julia_abcd.tar.gz

.PHONY: build
build:
	docker-compose build
