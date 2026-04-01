.PHONY: help cpp java clean

help:
	@echo "Targets:"
	@echo "  make cpp   Build the paired VSEARCH binary into build/cpp"
	@echo "  make java  Compile the Java launcher into build/java"
	@echo "  make clean Remove local build artifacts under build/"

cpp:
	./scripts/build_cpp.sh

java:
	./scripts/build_java.sh

clean:
	rm -rf build
