/*
 * MIT License
 *
 * Copyright (c) 2025 Betuel Sevindik, Maximilian Rottmann, Felix Baensch, Jonas
 * Schaub, Christoph Steinbeck, and Achim Zielesny
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.fragment.fingerprint;

import org.junit.jupiter.api.Test;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fragment.ExhaustiveFragmenter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Test class with usage examples for the fragment fingerprinter functionality.
 *
 * @version 1.1.0.0
 * @author Jonas Schaub, Maximilian Rottmann
 */
public class ExampleUsageTest {
    /**
     * Demonstrates the getFragmentsComponentsFloatMatrix() method with two fingerprinting modes:
     * bit array (binary presence/absence) and count array (fragment frequencies).
     *
     * The method processes an array of fragment-frequency maps and populates a pre-initialized
     * float matrix with fingerprint values. This is useful for batch processing multiple molecules
     * or generating matrix-based representations suitable for machine learning pipelines.
     *
     * The fragments used represent the 10 most frequently occurring functional group fragments from a COCONUT
     * database analysis of 1000 molecules.
     *
     * @see FragmentFingerprinter#getFragmentsComponentsFloatMatrix(Map[], float[][], boolean)
     * @see CountFingerprint for alternative fingerprinting approaches
     */
    @Test
    public void floatMatrixExampleUsageTest() {
        System.out.println("Float matrix example usage test:");
        //
        // SETUP:
        //
        // Master vector (pre-defined fingerprint) defines the standardized fragment set.
        // Each index in this list maps to a column position in the resulting fingerprint matrix.
        // Only fragments present in this master list will be included in the fingerprint.
        List<String> tmpFingerprintMasterList = new ArrayList<>(10);
        // SMILES representations of fragments (canonical/unique form required for consistent fingerprinting)
        // The order here determines the bit/column index mapping in the fingerprinter
        tmpFingerprintMasterList.add("C");
        tmpFingerprintMasterList.add("CC");
        tmpFingerprintMasterList.add("[H]OC");
        tmpFingerprintMasterList.add("*n(*)*");
        tmpFingerprintMasterList.add("*O*");
        tmpFingerprintMasterList.add("CCC");
        tmpFingerprintMasterList.add("C=C");
        tmpFingerprintMasterList.add("c");
        tmpFingerprintMasterList.add("*Cl");
        tmpFingerprintMasterList.add("CCCC");
        // Initialize FragmentFingerprinter with the master fragment list. This creates an internal
        // SMILES-to-position mapping (Map<String, Integer>) that will be used to:
        // 1) Determine fingerprint size (equal to number of unique fragments)
        // 2) Map fragment SMILES to their column indices in the output matrix
        // 3) Filter unrecognized fragments when generating fingerprints
        FragmentFingerprinter tmpFFp = new FragmentFingerprinter(tmpFingerprintMasterList);
        // Fragment-frequency map for a single molecule (represents one row in output matrix).
        // Key: fragment SMILES; Value: occurrence count in the molecule.
        // The getFragmentsComponentsFloatMatrix() method expects an array of such maps,
        // where each array element represents one molecule.
        Map<String, Integer> tmpFragmentsFrequenciesMap = new HashMap<>(16, 0.75f);
        // Populate with example fragment frequencies (sorted descending by abundance).
        // These represent the counts of each fragment found in a sample molecule.
        // Values will be converted to floats, 1.0f/0.0f in bit-array mode, or the retained values in count-array mode.
        tmpFragmentsFrequenciesMap.put("C", 10);
        tmpFragmentsFrequenciesMap.put("CC", 9);
        tmpFragmentsFrequenciesMap.put("[H]OC", 8);
        tmpFragmentsFrequenciesMap.put("*n(*)*", 7);
        tmpFragmentsFrequenciesMap.put("*O*", 6);
        tmpFragmentsFrequenciesMap.put("CCC", 5);
        tmpFragmentsFrequenciesMap.put("C=C", 4);
        tmpFragmentsFrequenciesMap.put("c", 3);
        tmpFragmentsFrequenciesMap.put("*Cl", 2);
        tmpFragmentsFrequenciesMap.put("CCCC", 1);
        // Additional fragments (*F, *I) are included to demonstrate the filtering mechanism:
        // getFragmentsComponentsFloatMatrix() only processes fragments that exist in the master list
        // (initialized in FragmentFingerprinter constructor). Unknown fragments are silently ignored.
        // This is why the output only shows values for the 10 master fragments, not these 2 additions.
        tmpFragmentsFrequenciesMap.put("*F", 1);
        tmpFragmentsFrequenciesMap.put("*I", 1);
        //
        // USAGE:
        //
        // Initialize two float matrices with different column counts:
        // - tmpDataMatrix: exactly matches fingerprint size (10 columns) - demonstrates standard usage
        // - tmpDataMatrixWithOverhang: larger than fingerprint size (15 columns) - tests matrix validation
        // The method validates that matrix columns >= fingerprint size; excess columns remain untouched.
        float[][] tmpDataMatrix = new float[1][10];
        float[][] tmpDataMatrixWithOverhang = new float[1][15];
        // Initialize array of fragment-frequency maps. In this example, a single-element array
        // is used (representing one molecule), but the method can process multiple molecules in batch.
        Map<String, Integer>[] tmpFragmentsMapsArray = new HashMap[1];
        tmpFragmentsMapsArray[0] = tmpFragmentsFrequenciesMap;
        // Call getFragmentsComponentsFloatMatrix() with anUseBitArrayStatement=true to generate a BIT-ARRAY fingerprint.
        // Internal method flow:
        // 1) Validates matrix dimensions (rows >= maps array length, columns >= fingerprint size)
        // 2) Iterates through each map in the array and calls private getFloatFingerprint()
        // 3) For recognized fragments: sets matrix[row][column] = 1.0f (presence indicator)
        // 4) For unrecognized fragments: ignores them (no value set)
        // 5) Uninitialized positions (fragments not in this molecule) remain 0.0f
        tmpFFp.getFragmentsComponentsFloatMatrix(tmpFragmentsMapsArray, tmpDataMatrix, true);
        // Visualize resulting matrix
        System.out.println("Matrix without overhang columns:");
        for (float[] tmpFloatArray : tmpDataMatrix) {
            System.out.println(Arrays.toString(tmpFloatArray));
        }
        /* Output:
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
         */
        //
        // Demonstrate matrix overhang handling with the same bit-array mode operation:
        // This matrix has 15 columns but the fingerprint only uses 10 (based on 10 master fragments).
        // The getFragmentsComponentsFloatMatrix() method:
        // 1) Only writes to indices 0-9 (columns within fingerprint size)
        // 2) Leaves columns 10-14 untouched (preserving any pre-existing values)
        // This allows reusing matrix memory or combining multiple fingerprints in one matrix.
        tmpFFp.getFragmentsComponentsFloatMatrix(
                tmpFragmentsMapsArray, tmpDataMatrixWithOverhang, true);
        // Visualize matrix
        System.out.println("Matrix WITH overhang columns:");
        for (float[] tmpFloatArray : tmpDataMatrixWithOverhang) {
            System.out.println(Arrays.toString(tmpFloatArray));
        }
        /* Output:
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
         */
        //
        // Switch to COUNT-ARRAY fingerprint mode (anUseBitArrayStatement=false).
        // Instead of binary 1.0f/0.0f for presence/absence, this mode stores actual fragment frequencies:
        // - Recognized fragments: matrix[row][column] = frequency value from the input map
        // - Unrecognized fragments: ignored (no entry)
        // - Missing fragments: remain 0.0f (can represent zero occurrences)
        tmpFFp.getFragmentsComponentsFloatMatrix(tmpFragmentsMapsArray, tmpDataMatrix, false);
        // Visualize matrix
        System.out.println("Matrix without overhang columns BUT count array behavior:");
        for (float[] tmpFloatArray : tmpDataMatrix) {
            System.out.println(Arrays.toString(tmpFloatArray));
        }
        /* Output:
        [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0]
         */
        // A matrix with overhang works analogous as in the example above and will not be shown explicitly.
    }
    /**
     * The intended use case of the fragment fingerprinter functionality is to encode the presence and absence of
     * substructures in a given molecule that result from a molecular fragmentation study, i.e. the algorithmic
     * extraction of specific substructures from input molecules. These substructures are automatically extracted and
     * can be represented by different string-based molecular structure encodings, like SMILES or InChI. Other
     * key-based substructure fingerprint functionalities require SMARTS strings as inputs and are therefore not
     * as ubiquitously applicable as the fragment fingerprint for this purpose.
     *
     * In the following, a molecular structure data set is imported that contains 100 natural products with a
     * naphthalene substructure taken from the COCONUT natural products database. These are fragmented using the CDK
     * ExhaustiveFragmenter functionality that breaks single non-ring bonds in input molecules to generate fragments.
     * The resulting fragments are collected together with their fraquencies as unique SMILES representations.
     * Fragments that occur more than two times are then used to initialise the fragment fingerprinter. At the end,
     * the "naphthalene-derivatives exhaustive fragmenter fingerprint" is generated for 3-hydroxy-2-naphthoic acid.
     */
    @Test
    public void chemicalExampleUsageTest() throws Exception {
        InputStream tmpInputStream = ExampleUsageTest.class.getResourceAsStream("coconut_naphthalene_substructure_search_result.sdf");
        //note: for the tutorial, make it InputStream tmpInputStream = new FileInputStream("\\path\\to\\coconut_naphthalene_substructure_search_result.sdf");
        IteratingSDFReader tmpSDFReader = new IteratingSDFReader(tmpInputStream, SilentChemObjectBuilder.getInstance());
        //This fragmentation scheme simply breaks single non-ring bonds.
        ExhaustiveFragmenter tmpFragmenter = new ExhaustiveFragmenter();
        //Default would be 6 which is too high for the short side chains in the input molecules
        tmpFragmenter.setMinimumFragmentSize(1);
        //ExhaustiveFragmenter has a convenience method .getFragments() that returns the generated fragments already as
        // unique SMILES strings, but to be explicit here, the fragments are retrieved as atom containers and unique
        // SMILES strings created in a second step. Also note that any other string-based molecular structure representation
        // like InChI could be used instead, but it should be canonical.
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        HashMap<String, Integer> tmpFrequenciesMap = new HashMap<>(50, 0.75f);
        while (tmpSDFReader.hasNext()) {
            IAtomContainer tmpMolecule = tmpSDFReader.next();
            tmpFragmenter.generateFragments(tmpMolecule);
            IAtomContainer[] tmpFragments = tmpFragmenter.getFragmentsAsContainers();
            for (IAtomContainer tmpFragment : tmpFragments) {
                String tmpSmilesCode = tmpSmiGen.create(tmpFragment);
                if (tmpFrequenciesMap.containsKey(tmpSmilesCode)) {
                    tmpFrequenciesMap.put(tmpSmilesCode, tmpFrequenciesMap.get(tmpSmilesCode) + 1);
                } else {
                    tmpFrequenciesMap.put(tmpSmilesCode, 1);
                }
            }
        }
        //Printing size of fragment set and all the fragment SMILES with their frequencies
        System.out.println(tmpFrequenciesMap.keySet().size());
        for (String tmpFragmentSmilesCode : tmpFrequenciesMap.keySet()) {
            System.out.println(tmpFragmentSmilesCode + ": " + tmpFrequenciesMap.get(tmpFragmentSmilesCode));
        }
        /*
         * Output:
         * 28
         * BrC1=CC=CC=2C=CC=CC12: 4
         * BrC=1C=CC2=CC(O)=CC=C2C1: 1
         * OC1=C[CH](OC)=CC=2C=CC=CC12: 1
         * BrC1=CC=CC2=[C]C=CC=C12: 1
         * BrC1=CC=CC=2C=[C]C=CC12: 1
         * O=CC: 1
         * O[NH](O)[CH]1=CC=CC=2C=CC=CC21: 1
         * O=CCl: 1
         * OC=1C=CC=2C=CC=CC2C1: 6
         * BrC1=CC=CC=2C=C(C=CC12)C: 1
         * ON=[CH3]: 2
         * ONO: 3
         * OC1=CC=CC=2C=CC=CC12: 4
         * NC1=CC=CC=2C=CC=CC12: 1
         * O=C[CH]1=CC=CC=2C=CC=CC21: 2
         * O=CO: 8
         * ON=C: 1
         * O=[S](=O)O: 5
         * C=1C=CC=2C=CC=CC2C1: 20
         * C=1C=CC=2C=C(C=CC2C1)C: 2
         * OC1=CC=CC=2C1=CC=CC2C: 1
         * O=N[CH]1=CC=C(O)C=2C=CC=CC21: 1
         * OC=1C=2C=CC=CC2C=CC1C: 1
         * [CH2][CH]=1C=CC=2C=CC=CC2C1: 1
         * BrC1=CC=CC=2C1=CC=CC2C: 1
         * OC1=CC=C(O)C=2C=CC=CC12: 2
         * C=1C=CC2=C(C1)C=CC=C2C: 1
         * O=COC: 1
         */
        //Collecting fragments that appear at least 2 times
        List<String> tmpFragmentsList = new ArrayList<>(28);
        for (String tmpFragment : tmpFrequenciesMap.keySet()) {
            if (tmpFrequenciesMap.get(tmpFragment) > 2) {
                tmpFragmentsList.add(tmpFragment);
            }
        }
        //Initialising fingerprinter
        FragmentFingerprinter tmpNaphthaleneFingerprinter = new FragmentFingerprinter(tmpFragmentsList);
        System.out.println(tmpNaphthaleneFingerprinter.getSize());
        /*
         * Output: 7
         *
         * Only 7 out of the 28 fragments appear more than 2 times and are included in the fingerprint (see above).
         */
        //Parsing 3-hydroxy-2-naphthoic acid, fragmenting it, and creating its fingerprint
        String tmpCNP0437667SmilesString = "O=C(O)C1=CC=2C=CC=CC2C=C1O"; //3-hydroxy-2-naphthoic acid
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpFragmenter.generateFragments(tmpSmiPar.parseSmiles(tmpCNP0437667SmilesString));
        IAtomContainer[] tmpFragments = tmpFragmenter.getFragmentsAsContainers();
        List<String> tmpCNP0437667Fragments = new ArrayList(10);
        for (IAtomContainer tmpFragment : tmpFragments) {
            tmpCNP0437667Fragments.add(tmpSmiGen.create(tmpFragment));
        }
        IBitFingerprint tmpCNP0437667BitFP = tmpNaphthaleneFingerprinter.getBitFingerprint(tmpCNP0437667Fragments);
        for (int i = 0; i < tmpNaphthaleneFingerprinter.getSize(); i++) {
            System.out.println(tmpNaphthaleneFingerprinter.getBitDefinition(i) + ": " + tmpCNP0437667BitFP.get(i));
        }
        /*
         * Output:
         * BrC1=CC=CC=2C=CC=CC12: false
         * OC=1C=CC=2C=CC=CC2C1: true
         * ONO: false
         * OC1=CC=CC=2C=CC=CC12: false
         * O=CO: true
         * O=[S](=O)O: false
         * C=1C=CC=2C=CC=CC2C1: false
         *
         * 3-hydroxy-2-naphthoic acid contains the formic acid and the naphthol fragments. It does not produce a
         * naphthalene fragment because the hydroxy fragment is too small to be considered on its own, according to the CDK
         * ExhaustiveFragmenter.
         */
    }
}
