.PHONY: help cpp java stage-install clean

help:
	@echo "Targets:"
	@echo "  make cpp   Build the paired VSEARCH binary into build/cpp"
	@echo "  make java  Compile the Java launcher into build/java"
	@echo "  make stage-install PREFIX=/path  Stage installable files into a prefix"
	@echo "  make clean Remove local build artifacts under build/"

cpp:
	./scripts/build_cpp.sh

java:
	./scripts/build_java.sh

stage-install:
	./scripts/stage_install.sh "$(PREFIX)"

clean:
	rm -rf build
