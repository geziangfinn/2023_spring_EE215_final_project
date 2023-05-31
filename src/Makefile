core_router: core_router.cpp# Build the core_router
	g++ -o core_router ../src/core_router.cpp

run: # Run the core router for all benchmarks
	./core_router ../benchmarks bench1 ../out > ../out/bench1
	./core_router ../benchmarks bench2 ../out > ../out/bench2
	./core_router ../benchmarks bench3 ../out > ../out/bench3
	./core_router ../benchmarks bench4 ../out > ../out/bench4
	./core_router ../benchmarks bench5 ../out > ../out/bench5
	./core_router ../benchmarks fract2 ../out > ../out/fract2
	./core_router ../benchmarks primary1 ../out > ../out/primary1
	./core_router ../benchmarks industry1 ../out > ../out/industry1

.PHONY:
clean: # Clean the build directory, and routing results
	rm -f core_router
	rm -f ../out/*.route