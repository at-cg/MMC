MMC
=
MMC is a open-source program for disk-based counting of (w,k)-minimizers from (possibly gzipped) FASTQ/FASTA files. Given parameters k and w (default w = k), MMC samples the minimizers and give the statistics of their count in an efficient manner. Many applications (such as, GWAS, BC-dissimilarity computation, de novo Genome Size Estimation ...) that requires the count of all the k-mers in a given sequencing file, can be performed using a fraction of time and memory using only the minimizer statistics, which MMC provides the user. This program has been built on top of [KMC3](https://github.com/refresh-bio/KMC).

Benchmark
=
![image](https://user-images.githubusercontent.com/40889593/141755297-d2ff9e99-7d01-4c3c-9541-fbd316e2168c.png)
We have benchmarked MMC with two of the most common tools for counting all k-mers. We can clearly observe from the above plots that counting (w,k)-minimizers using MMC is around 2x faster than KMC3 and 10x faster than Jelllyfish2 (JF2). MMC is also 2x space efficient compared to KMC3. But, the best part about using MMC is that is only samples minimizers, which results in subsequent downstream operations an order of magnitude faster in minimizer-space compared to the classical all k-mer counterparts.
Installation
=
```sh
git clone https://github.com/at-cg/MMC.git
cd MMC
```
After that, you can run ```make``` to compile mmc, mmc_dump, and mmc_tools applications.

**Note:** Some parts of MMC use C++14 features, so you need a compatible C++ compiler, e.g., gcc 4.9+ or clang 3.4+

Binaries
=
After compilation you will obtain two binaries:
* bin/mmc - the main program for counting minimizer occurrences
* bin/mmc_dump - the program listing minimizers in a database produced by mmc
* bin/mmc_tools - the program allowing to manipulate mmc databases (set operations, transformations, etc.)

Script
=
The scripts folder contains:
* bc_script.cpp - the script to compute Bray-Curtis dissimilarity between two metagenomes.

Usage: 
```sh
g++ bcscript.cpp -o main && ./main <input1> <input2>
```


License
=
* MMC software distributed under GNU GPL 3 licence.

* libbzip2 is open-source (BSD-style license)

* gzip is free, open-source

* pybind11 (https://github.com/pybind/pybind11) is open-source (BDS-style license)

In case of doubt, please consult the original documentations.

Warranty
=
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING
THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

