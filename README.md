MMC
=
MMC is a disk-based programm for counting (w,k)-minimizers from (possibly gzipped) FASTQ/FASTA files.

Installation
=
```sh
git clone https://github.com/at-cg/MMC.git
cd MMC
```
If needed, you can also redefine maximal length of k-mer, which is 256 in the current version.

Note: MMC is highly optimized and spends only as many bytes for k-mer (rounded up to 8) as
necessary, so using large values of MAX_K does not affect the KMC performance for short k-mers.

Some parts of MMC use C++14 features, so you need a compatible C++ compiler, e.g., gcc 4.9+ or clang 3.4+

After that, you can run ```make``` to compile mmc and mmc_dump applications.

Directory structure
=
 * bin           - main directory of MMC (programs after compilation will be stored here)
 * kmer_counter  - source code of kmc program
 * kmer_counter/libs - compiled binary versions of libraries used by KMC
 * kmc_api       - C++ source codes implementing API; must be used by any program that wants to process databases produced by kmc
 * kmc_dump      - source codes of kmc_dump program listing k-mers in databases produced by kmc
 * py_kmc_api    - python wrapper for kmc API

Binaries
=
After compilation you will obtain two binaries:
* bin/mmc - the main program for counting k-mer occurrences
* bin/mmc_dump - the program listing k-mers in a database produced by mmc
* bin/mmc_tools - the program allowing to manipulate mmc databases (set operations, transformations, etc.)


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

