// main test cfsem-rs built with cpp

#include <cstdio>

// include rust function signatures
// rust functions need
// #[no_mangle]
// and pub extern "C" fn
// signatures to be callable
extern "C" {
double ellipk(double);
}

int main() {

	printf("Initialize cfsem-rs\n");
	printf("cfsem ellipk() function: (%0.3f)\n", ellipk(0.1));

	// success
	return 0;
}
