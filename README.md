# F2-Cmax-revisited

## Description
This project implements the Johnson's rule to solve F2-Cmax scheduling problem as well as the revisited version described in https://doi.org/10.48550/arXiv.2512.06119 

## Installation

### Dependencies

* Download and install **nlohmann/json** version 3.11.3 from GitHub: https://github.com/nlohmann/json/archive/refs/tags/v3.11.3.zip
  + Place the downloaded archive in `./lib/` directory.
  + Unzip the archive using the command `unzip v3.11.3.zip`.

The project uses CMake for installation. To generate and build the project, follow these steps:

* Navigate to the project root directory.
* Run the following command to generate the project in release mode: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`
* Run the following command to generate the project in debug mode: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug`
* Change into the `/build` directory and run the command `make`.

The compiled binary program will be located in either `/bin/release` or `/bin/debug`, depending on the build mode chosen.

**Note:** The project was not compiled on Windows, some error may appear. Please feel free to open an issue if you encounter any issues during the compilation process.

## Specification

### Instances

An instance file must be in "txt" format and follow the specific pattern:
```
name:Instance Name
n:Integer # Number of jobs
Jobs:
...
pj1 \t pj2 
...
```
In the job description, the processing time `pj1` is the processing on machine M1, and `pj2` is the one on machine M2 must be separated by a `\t` character.

## Usage

To run the F2-Cmax-revisited project, execute the binary program located in either `./bin/release` or `./bin/debug` directory. You must provide a configuration file. For reference, you can find two example files: `example_config_generate.json` and `example_config_solve.json`.

All instances used in the paper are in "instances.tar.bz2"

### Configuration File (JSON)

The project uses JSON files as configuration files. Below, we describe the structure of these files and how to use them.

#### Generate

To generate instances, you'll need to provide a JSON configuration file with the following structure:
```
{
    "generate": {
        // Seed for generating instances
        "seed": <int>,
        
        // Parameters for generating one instance
        "instances": [
            {
                // Base path where to save instances. It creates the directory if it does not exist.
                "basePath": "<string>",
                
                // Number of instances to generate
                "numberInstance": <int>,
                
                "paramInstance": {
                    // The number of jobs that the leader has to select.
                    "n": <int>,
                    
                    // Total number of jobs in the instance
                    "distribution": <string> (uniform,negative_binomial,geometric,poisson,exponential,gamma,weibull,lognormal),
                    
                    // Processing time range for each job
                    "pi": {
                        "inf": <int>, // Lower bound of processing time
                        "sup": <int>  // Upper bound of processing time
                    },
                }
            }
        ]
    }
}
```
This configuration file specifies the parameters needed to generate instances, including the seed value, instance base path, and parameters for generating a single instance. The `paramInstance` object contains additional parameters that control the generation process, such as the number of jobs, and the ranges for processing times, the distribution to follow.

#### Solve

To solve instances, there are different parameters described below:

```
"solve": {
    // Level of verbose
    "verbose": <int>,
    // Path where to save the results. If directory does not exist, it will be created.
    "output": "<string>",
    // List of methods used to define parameters for each method. This is described below.
    "methods": [
        {
            // Level of verbose
            "verbose": <int>,
            // use the revisited algorithm
            "useRevisited": <bool>, 
            // List of instances to solve. Each object is composed of only one attribute:
            "instances": [
                {
                    // Path to the instance file to solve
                    "path": "<string>"
                }
            ]
        }
    ]
}
```

## Contributing

The main contributor is Quentin SCHAU. If you want to contribute to this project, you should reach out to Quentin SCHAU at quentin.schau@univ-tours.fr or quentin.schau@polito.it .

## Authors and acknowledgment

This code implements the problem and the revisited algorithm defined in  	
https://doi.org/10.48550/arXiv.2512.06119

## License
This project is under GPU Licence. 

