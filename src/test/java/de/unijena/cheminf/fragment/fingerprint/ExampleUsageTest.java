/*
 * MIT License
 *
 * Copyright (c) 2023 Betuel Sevindik, Felix Baensch, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
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
 *
 */

package de.unijena.cheminf.fragment.fingerprint;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
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
import java.util.HashMap;
import java.util.List;

/**
 * Test class with usage examples for the fragment fingerprinter functionality.
 *
 * @version 1.0.0.0
 * @author Jonas Schaub
 */
public class ExampleUsageTest {
    //<editor-fold desc="private static class variables" defaultstate="collapsed">
    /**
     * List of all key fragments passed during fingerprint initialization.
     */
    private static List<String> fragmentsList;
    /**
     * FragmentFingerprinter
     */
    private static FragmentFingerprinter naphthaleneFingerprinter;
    /**
     * molecule fragments of naphthalene derivate
     */
    private static List<String> cNP0437667Fragments;
    //</editor-fold>
    //
    //<editor-fold desc="BeforeAll method" defaultstate="collapsed">
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
     * The resulting fragments are collected together with their frequencies as unique SMILES representations.
     * Fragments that occur more than two times are then used to initialise the fragment fingerprinter.
     *
     */
    @BeforeAll
    public static void SetUpChemicalExampleUsageTest() throws CDKException {
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
        System.out.println("Printing size of the fragment set and all the fragment SMILES with their frequencies");
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
        ExampleUsageTest.fragmentsList = new ArrayList<>(28);
        for (String tmpFragment : tmpFrequenciesMap.keySet()) {
            if (tmpFrequenciesMap.get(tmpFragment) > 2) {
                ExampleUsageTest.fragmentsList.add(tmpFragment);
            }
        }
        //Parsing 3-hydroxy-2-naphthoic acid, fragmenting it, and creating its fingerprint
        String tmpCNP0437667SmilesString = "O=C(O)C1=CC=2C=CC=CC2C=C1O"; //3-hydroxy-2-naphthoic acid
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpFragmenter.generateFragments(tmpSmiPar.parseSmiles(tmpCNP0437667SmilesString));
        IAtomContainer[] tmpFragments = tmpFragmenter.getFragmentsAsContainers();
        ExampleUsageTest.cNP0437667Fragments = new ArrayList(10);
        for (IAtomContainer tmpFragment : tmpFragments) {
            ExampleUsageTest.cNP0437667Fragments.add(tmpSmiGen.create(tmpFragment));
        }
    }
    //
    /**
     * The "naphthalene-derivatives exhaustive fragmenter fingerprint" is generated for 3-hydroxy-2-naphthoic acid.
     * The size of the generated fingerprint is checked.
     */
    @Test
    public void chemicalExampleUsageTest() throws Exception {
        //Initialising fingerprinter and testing the fingerprint size
        FragmentFingerprinter tmpNaphthaleneFingerprinter = new FragmentFingerprinter(ExampleUsageTest.fragmentsList);
        int tmpFingerprintSizeOfNaphthalene = 7;
        /*
         * Output: 7 is
         *
         * Only 7 out of the 28 fragments appear more than 2 times and are included in the fingerprint (see above).
         */
        Assertions.assertEquals(tmpFingerprintSizeOfNaphthalene, tmpNaphthaleneFingerprinter.getSize());
    }
    //
    /**
     * The bit definition in position 5 in the naphthalene fingerprint is checked.
     */
    @Test
    public void chemicalExampleUsageTestBitDefinitionInNaphthaleneFingerprint() {
        FragmentFingerprinter tmpNaphthaleneFingerprinter = new FragmentFingerprinter(ExampleUsageTest.fragmentsList);
        // In position 5 in the fingerprint the fragment O=[S](=O)O is present
        String tmpBitDefinitionInFingerprint = "O=[S](=O)O";
        Assertions.assertEquals(tmpBitDefinitionInFingerprint, tmpNaphthaleneFingerprinter.getBitDefinition(5));
    }
    //
    /**
     * The method generates the bit fingerprint for the given molecule.
     * In the method is not tested. The method is intended only as a small illustration of the naphthalene fingerprint.
     */
    @Test
    public void illustrateSomeResultsOfChemicalExampleUsage() {
        FragmentFingerprinter tmpNaphthaleneFingerprinter = new FragmentFingerprinter(ExampleUsageTest.fragmentsList);
        IBitFingerprint tmpCNP0437667BitFP = tmpNaphthaleneFingerprinter.getBitFingerprint(ExampleUsageTest.cNP0437667Fragments);
        System.out.println("\n");
        System.out.println("Print of the bit definition in the fingerprint and the associated value. The value is true if the bit with the given index is currently set in this fingerprint; otherwise, the result is false.");
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
    /**
     * At the very basic level, the fragment fingerprinter is a simple string matching-based functionality for creating
     * bit and count vectors based on a set of initialisation strings that the fingerprinter checks given sets of
     * strings for.
     *
     * The following illustrates this with a basic example.
     */
    @Test
    public void generalExampleUsageTest() throws Exception {
        List<String> tmpInitialisationStrings = List.of("Hannah", "Sam", "John", "Hugo", "Tim");
        //Initialising the fingerprinter with the list of names
        FragmentFingerprinter tmpFingerprinter = new FragmentFingerprinter(tmpInitialisationStrings);
        //Creating a set of names to generate a fingerprint for
        List<String> tmpMyPartyPeople = List.of("Hugo", "Hannah", "Sam", "Maria");
        //Generating the bit fingerprint
        IBitFingerprint tmpMyPartyPeopleFP = tmpFingerprinter.getBitFingerprint(tmpMyPartyPeople);
        int tmpCardinalityCheck = 3;
        //execution of the test; is checked whether the number of positive bits is correct.
        Assertions.assertEquals(tmpCardinalityCheck, tmpMyPartyPeopleFP.cardinality());
        // In the following lines, some results are printed for illustration purposes:
        System.out.println("\n");
        System.out.println("General example usage of the fragment fingerprinter. The following illustrates this with a basic example.");
        System.out.println();
        System.out.println(tmpMyPartyPeopleFP.cardinality());
        /*
         * Output: 3
         *
         * I.e. 3 positive bits in the fingerprint. Maria is ignored since she was not part of the initialisation set.
         */
        System.out.println(tmpMyPartyPeopleFP.asBitSet().toString());
        /*
         * Output: {0, 1, 3}
         *
         * Hannah is represented by position 0, Sam by position 1, and Hugo by position 3 in
         * the fingerprint and these positions are positive in the bit fingerprint.
         */
        //Printing the bit definitions and whether they are positive or not in the "party" set of names
        for (int i = 0; i < tmpFingerprinter.getSize(); i++) {
            System.out.println(tmpFingerprinter.getBitDefinition(i) + ": " + tmpMyPartyPeopleFP.get(i));
        }
        /*
         * Output:
         * Hannah: true
         * Sam: true
         * John: false
         * Hugo: true
         * Tim: false
         */
    }
}
