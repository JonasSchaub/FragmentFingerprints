[![DOI](https://zenodo.org/badge/548827010.svg)](https://zenodo.org/doi/10.5281/zenodo.7998215)
[![Javadoc](https://img.shields.io/badge/JavaDoc-Online-green)](https://jonasschaub.github.io/FragmentFingerprints/javadoc/latest/index.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-blue.svg)](https://GitHub.com/JonasSchaub/FragmentFingerprints/graphs/commit-activity)
[![build](https://github.com/JonasSchaub/FragmentFingerprints/actions/workflows/gradle.yml/badge.svg)](https://github.com/JonasSchaub/FragmentFingerprints/actions/workflows/gradle.yml)
[![GitHub issues](https://img.shields.io/github/issues/JonasSchaub/FragmentFingerprints.svg)](https://GitHub.com/JonasSchaub/FragmentFingerprints/issues/)
[![GitHub contributors](https://img.shields.io/github/contributors/JonasSchaub/FragmentFingerprints.svg)](https://GitHub.com/JonasSchaub/FragmentFingerprints/graphs/contributors/)
[![GitHub release](https://img.shields.io/github/release/JonasSchaub/FragmentFingerprints.svg)](https://github.com/JonasSchaub/FragmentFingerprints/releases/)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/io.github.jonasschaub/FragmentFingerprints/badge.svg)](https://maven-badges.herokuapp.com/maven-central/io.github.jonasschaub/FragmentFingerprints)
[![Quality Gate Status](https://sonarcloud.io/api/project_badges/measure?project=JonasSchaub_FragmentFingerprints&metric=alert_status)](https://sonarcloud.io/summary/new_code?id=JonasSchaub_FragmentFingerprints)
# FragmentFingerprints
A library to generate fingerprints for molecular structures based on a set of fragments

## Description
The library generates fragment fingerprints based on pre-defined fragments, 
which can be set by the user, and can generate both bit and count fingerprints. 
Fragment fingerprints are created by matching fragments or substructures of a 
molecule with pre-defined fragments. If a match is found, the corresponding positions 
in the fingerprint are filled. The special feature of the fragment fingerprinter is that fingerprints 
are generated exclusively by comparing unique SMILES (Strings). This means that both the pre-defined fragments and 
the substructures or fragments of the molecule for which the fingerprint is 
being generated must be represented as unique SMILES strings. The implementation of the fragment fingerprinter is based on
the Chemistry Development Kit (CDK).

## Contents of this repository
### Sources
The <a href="https://github.com/JonasSchaub/FragmentFingerprints/tree/main/src">"src"</a> subfolder contains
all source code packages including JUnit tests.

### Tests
The test class <i>FragmentFingerprinterTest</i> tests the functionalities of fragment fingerprinter.
Among other things, it tests whether the bit and count fingerprint of a molecule has been generated 
correctly. Furthermore, various methods of the CountFingerprint and 
BitSetFingerprint classes are tested.

### Test resources
The test <a href="https://github.com/JonasSchaub/FragmentFingerprints/tree/FragmentFingerprint/src/test/resources/de/unijena/cheminf/fragment/fingerprint">"resources"</a> subfolder
contains two text files. The text file named "FragmentList.txt" contains all key fragments. And the file 
named "MoleculeList.txt" contains fragments/substructures of molecules. 
In total, 10 molecules with their corresponding fragments are stored in the file.

### Performance Test CMD Application
The folder <a href="https://github.com/JonasSchaub/FragmentFingerprints/tree/FragmentFingerprint/PerformanceTestCMDApplication">"PerformanceTestCMDApplication"</a>
contains the executable JAVA archive <i>FragmentFingerprints-fat.jar</i>.
It can be executed from the command-line (command: java -jar) to do a performance snapshot of fragment fingerprinter's scaling behaviour for
a growing number of input molecules.
For more details see the file <a href="https://github.com/JonasSchaub/FragmentFingerprints/blob/FragmentFingerprint/PerformanceTestCMDApplication/Performance_test_instruction.txt">"Performance_test_instruction.txt"</a>

## Example initialization and usage of the FragmentFingerprinter
see in <a href="https://github.com/JonasSchaub/FragmentFingerprints/wiki">"wiki"</a>

## Installation
FragmentFingerprints is hosted as a package/artifact on the sonatype maven central repository. See the
<a href="https://central.sonatype.com/artifact/io.github.jonasschaub/FragmentFingerprints/">artifact page</a> for installation
guidelines using build tools like maven or gradle.
<br>
To install FragmentFingerprints via its JAR archive, you can get it from the
<a href="https://github.com/JonasSchaub/FragmentFingerprints/releases">releases</a>. Note that other dependencies
will need to be installed via JAR archives as well this way.
<br>
In order to open the project locally, e.g. to extend it, download or clone the repository and
open it in a Gradle-supporting IDE (e.g. IntelliJ) as a Gradle project and execute the build.gradle file.
Gradle will then take care of installing all dependencies. A Java Development Kit (JDK) of version 17 or higher must also
be pre-installed.

## Dependencies
**Needs to be pre-installed:**
* Java Development Kit (JDK) version 17
  * [Adoptium OpenJDK](https://adoptium.net) (as one possible source of the JDK)
* Gradle version 7.3
  * [Gradle Build Tool](https://gradle.org)

**Managed by Gradle:**
* Chemistry Development Kit (CDK) version 2.9
  * [Chemistry Development Kit on GitHub](https://cdk.github.io/)
  * License: GNU Lesser General Public License 2.1
* JUnit Jupiter version 5.9.1
  * [JUnit ](https://junit.org/junit5/)
  * License: Eclipse Public License - v 2.0

## References and useful links
**Chemistry Development Kit (CDK)**
* [Chemistry Development Kit on GitHub](https://cdk.github.io/)
* [Steinbeck C, Han Y, Kuhn S, Horlacher O, Luttmann E, Willighagen EL. The Chemistry Development Kit (CDK): An Open-Source Java Library for Chemo- and Bioinformatics. J Chem Inform Comput Sci. 2003;43(2):493-500.](https://dx.doi.org/10.1021%2Fci025584y)
* [Steinbeck C, Hoppe C, Kuhn S, Floris M, Guha R, Willighagen EL. Recent Developments of the Chemistry Development Kit (CDK) - An Open-Source Java Library for Chemo- and Bioinformatics. Curr Pharm Des. 2006; 12(17):2111-2120.](https://doi.org/10.2174/138161206777585274)
* [May JW and Steinbeck C. Efficient ring perception for the Chemistry Development Kit. J. Cheminform. 2014; 6:3.](https://dx.doi.org/10.1186%2F1758-2946-6-3)
* [Willighagen EL, Mayfield JW, Alvarsson J, Berg A, Carlsson L, Jeliazkova N, Kuhn S, Pluska T, Rojas-Chertó M, Spjuth O, Torrance G, Evelo CT, Guha R, Steinbeck C, The Chemistry Development Kit (CDK) v2.0: atom typing, depiction, molecular formulas, and substructure searching. J Cheminform. 2017; 9:33.](https://doi.org/10.1186/s13321-017-0220-4)
* [Groovy Cheminformatics with the Chemistry Development Kit](https://github.com/egonw/cdkbook)