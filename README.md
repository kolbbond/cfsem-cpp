Leverage the [Corrosion](https://github.com/corrosion-rs/corrosion)
to build the cfsem-rs library with cmake

## Letting Rust functions be called from C
Add to your interface functions 
##
    #[no_mangle]
    pub extern "C" fn 

# Install and build
Download the repo with the cfsem-rs library 

##
    `git clone --recursive git@github.com:kolbbond/cfsem-cpp.git`

Build using the regular cmake process
##
    mkdir build;
    cd build;
    cmake ..;
    make;
