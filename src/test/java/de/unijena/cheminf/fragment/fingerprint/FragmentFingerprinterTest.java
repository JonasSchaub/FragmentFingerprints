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
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.fragment.ExhaustiveFragmenter;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;

/**
 * Class to test the correct working of FragmentFingerprinter
 *
 * @author Betuel Sevindik
 * @version 1.0.0.0
 */
public class FragmentFingerprinterTest {
    //<editor-fold desc="private static final class variables" defaultstate="collapsed">
    /**
     * Name of file for writing fingerprints results.
     */
    private static final String FINGERPRINTS_FILE_NAME = "Fingerprints";
    /**
     * Initial capacity value for molecule list
     */
    private static final int INITIAL_CAPACITY_VALUE = 10;
    /**
     *  Initial capacity value of the list that stores the key fragments.
     */
    private static final int INITIAL_CAPACITY_VALUE_FOR_FRAGMENT_LIST = 28;
    /**
     * Initial capacity value of maps.
     */
    private static final double INITIAL_CAPACITY_VALUE_FOR_MAP = 1.5;
    //</editor-fold>
    //
    //<editor-fold desc="private static class variables" defaultstate="collapsed">
    /**
     * The list contains a collection of fragments and their frequencies that were read from a
     * CSV file and stored as a HashMap.
     */
    private static ArrayList<HashMap<String, Integer>> moleculeFragmentList;
    /**
     * The list includes all key fragments that are set during the initialization of the fingerprinter.
     */
    private static ArrayList<String> fragmentList;
    /**
     * fragment fingerprinter
     */
    private static FragmentFingerprinter fragmentFingerprinter;
    /**
     * fragment fingerprinter
     */
    private static FragmentFingerprinter naphthaleneFingerprinter;
    /**
     * Bit fingerprint
     */
    private static IBitFingerprint bitFingerprintTest;
    /**
     * Bit fingerprint of naphthalene derivate.
     */
    private static IBitFingerprint cNP0437667BitFP;
    /**
     * Count fingerprint
     */
    private static CountFingerprint countFingerprintTest;
    /**
     * Fragments of Naphthalene dervivate.
     */
    private static List<String> cNP0437667Fragments;
    /**
     * List contains molecule fragments without duplicates
     */
   private static ArrayList<String> dataForGeneratingBitFingerprint;
    /**
     * List contains molecule fragments with desired fragment duplicates.
     */
   private static ArrayList<String> countListOfUniqueSmiles;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Empty Constructor
     */
    public FragmentFingerprinterTest() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="BeforeAll method" defaultstate="collapsed">
    /**
     * Start generating fingerprinting data.
     * To create the fingerprints, 2 text files are used here. One of these text files contains the predefined
     * fragments and the other text file contains the fragments associated with the molecules.
     * Structure of the text files can be seen in resources folder.
     *
     * "at"BeforeAll ensures that the setUp method is only executed once.
     *
     * @throws Exception is thrown if anything goes wrong.
     */
    @BeforeAll
    public static void setUp() throws Exception {
        BufferedReader tmpFragmentSetReader;
        BufferedReader tmpMoleculeFragmentsReader;
        tmpFragmentSetReader = new BufferedReader(new FileReader("src/test/resources/de/unijena/cheminf/fragment/fingerprint/FragmentList.txt"));
        tmpMoleculeFragmentsReader = new BufferedReader(new FileReader("src/test/resources/de/unijena/cheminf/fragment/fingerprint/MoleculeFragments.txt"));
        /* Read CSV file to obtain fragments used to create the fingerprint. A fragment in the form of
        unique SMILES is displayed at the first position of each line. Only these unique
        SMILES from the file are stored in the list.
         */
        FragmentFingerprinterTest.fragmentList = new ArrayList<>(FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE_FOR_FRAGMENT_LIST);
        FragmentFingerprinterTest.moleculeFragmentList = new ArrayList<>(FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE);
        String tmpLine;
        String tmpSeparatorComma = ",";
        while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
            String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
            FragmentFingerprinterTest.fragmentList.add(tmpSmilesOfFragments[0]);
        }
        tmpFragmentSetReader.close();
        // removing header line value
        FragmentFingerprinterTest.fragmentList.remove(0);
        // Read CSV file
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpListOfMoleculesFragmentsAndFrequenciesList = new ArrayList<>(FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE);
        String tmpMoleculeLine;
        while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
            String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
            tmpListOfMoleculesFragmentsAndFrequenciesList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
        }
        tmpMoleculeFragmentsReader.close();
        List<String> tmpNameAndMoleculeFragmentsAndFrequenciesList;
        FragmentFingerprinterTest.fragmentFingerprinter = new FragmentFingerprinter(FragmentFingerprinterTest.fragmentList);
        String tmpFingerprintOutputPath = (new File("").getAbsoluteFile().getAbsolutePath()) + File.separator;
        new File(tmpFingerprintOutputPath + FINGERPRINTS_FILE_NAME).mkdirs();
        File tmpFingerprintResultFile = new File(tmpFingerprintOutputPath + FINGERPRINTS_FILE_NAME +File.separator+ FragmentFingerprinterTest.FINGERPRINTS_FILE_NAME  + ".txt");
        FileWriter tmpFingerprintResultsFileWriter = new FileWriter(tmpFingerprintResultFile, false);
        PrintWriter tmpFingerprintResultPrintWriter = new PrintWriter(tmpFingerprintResultsFileWriter);
        int[] tmpBitArray;
        for (int i= 1; i < tmpListOfMoleculesFragmentsAndFrequenciesList.size(); i++) {
            tmpNameAndMoleculeFragmentsAndFrequenciesList = tmpListOfMoleculesFragmentsAndFrequenciesList.get(i);
            List<String> tmpListWithoutNameAndMoleculeSmiles = tmpNameAndMoleculeFragmentsAndFrequenciesList.subList(2, tmpNameAndMoleculeFragmentsAndFrequenciesList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>((int) (tmpListWithoutNameAndMoleculeSmiles.size()*FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE_FOR_MAP));
            FragmentFingerprinterTest.dataForGeneratingBitFingerprint = new ArrayList<>(FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE);
            for (int j = 0; j < tmpListWithoutNameAndMoleculeSmiles.size(); j++) {
                if (j % 2 == 0) { // magic number to store the fragment SMILES and their frequencies from the file into a HashMap
                    tmpMoleculeFragmentsMap.put(tmpListWithoutNameAndMoleculeSmiles.get(j), Integer.valueOf(tmpListWithoutNameAndMoleculeSmiles.get(j + 1)));
                    FragmentFingerprinterTest.dataForGeneratingBitFingerprint.add(tmpListWithoutNameAndMoleculeSmiles.get(j));
                }
            }
            //Illustration of the results of the bit arrays for the specified molecules
            tmpBitArray = FragmentFingerprinterTest.fragmentFingerprinter.getBitArray(FragmentFingerprinterTest.dataForGeneratingBitFingerprint);
            // number of characters
            int tmpLength = java.util.Arrays.toString(tmpBitArray).length();
            tmpFingerprintResultPrintWriter.println(java.util.Arrays.toString(tmpBitArray).substring(1, tmpLength - 1));
            tmpFingerprintResultPrintWriter.flush();
            // add all molecule maps
            FragmentFingerprinterTest.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
        }
        tmpFingerprintResultPrintWriter.close();
        // Objects necessary for the test are created (used only in @Test)
        FragmentFingerprinterTest.countFingerprintTest = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.moleculeFragmentList.get(FragmentFingerprinterTest.moleculeFragmentList.size() - 1));
        FragmentFingerprinterTest.bitFingerprintTest = FragmentFingerprinterTest.fragmentFingerprinter.getBitFingerprint(FragmentFingerprinterTest.dataForGeneratingBitFingerprint);

        // Variamycin fragments
        FragmentFingerprinterTest.countListOfUniqueSmiles = new ArrayList<>(35);
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*C(*)=O");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*C(*)=O");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]Oc");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]Oc");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*O*");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("c1cc(cc2ccc(cc12)CC(C)C)C");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("C");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("*O*");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("C");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("CCC");
        FragmentFingerprinterTest.countListOfUniqueSmiles.add("[H]OC");
        /*
         * In the following: Creation of a test for chemical use of the fragment fingerprinter.
         *
         * In the following, a molecular structure data set is imported that contains 100 natural products with a
         * naphthalene substructure taken from the COCONUT natural products database. These are fragmented using the CDK
         * ExhaustiveFragmenter functionality that breaks single non-ring bonds in input molecules to generate fragments.
         * The resulting fragments are collected together with their fraquencies as unique SMILES representations.
         * Fragments that occur more than two times are then used to initialise the fragment fingerprinter. At the end,
         * the "naphthalene-derivatives exhaustive fragmenter fingerprint" is generated for 3-hydroxy-2-naphthoic acid.
         */
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
        //Collecting fragments that appear at least 2 times
        List<String> tmpFragmentsList = new ArrayList<>(28);
        for (String tmpFragment : tmpFrequenciesMap.keySet()) {
            if (tmpFrequenciesMap.get(tmpFragment) > 2) {
                tmpFragmentsList.add(tmpFragment);
            }
        }
        // Initialising fingerprinter
        FragmentFingerprinterTest.naphthaleneFingerprinter = new FragmentFingerprinter(tmpFragmentsList);
        //Parsing 3-hydroxy-2-naphthoic acid, fragmenting it, and creating its fingerprint
        String tmpCNP0437667SmilesString = "O=C(O)C1=CC=2C=CC=CC2C=C1O"; //3-hydroxy-2-naphthoic acid
        SmilesParser tmpSmiPar = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpFragmenter.generateFragments(tmpSmiPar.parseSmiles(tmpCNP0437667SmilesString));
        IAtomContainer[] tmpFragments = tmpFragmenter.getFragmentsAsContainers();
        FragmentFingerprinterTest.cNP0437667Fragments = new ArrayList(10);
        for (IAtomContainer tmpFragment : tmpFragments) {
            FragmentFingerprinterTest.cNP0437667Fragments.add(tmpSmiGen.create(tmpFragment));
        }
        // create bit fingerprint
        FragmentFingerprinterTest.cNP0437667BitFP = FragmentFingerprinterTest.naphthaleneFingerprinter.getBitFingerprint(FragmentFingerprinterTest.cNP0437667Fragments);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Test methods" defaultstate="collapsed">
    /**
     * Tests the number of positive bits in the bit fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void cardinalityTest() {
        int tmpNumberPositiveBitsTest = 9;
        int tmpNumberOfPositiveBitsVariamycin = FragmentFingerprinterTest.bitFingerprintTest.cardinality();
        Assertions.assertEquals(tmpNumberPositiveBitsTest, tmpNumberOfPositiveBitsVariamycin);
    }
    //
    /**
     * Tests the size of a bit fingerprint.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void bitFingerprintSizeTest() {
        long tmpVariamycinFingerprintSize =  FragmentFingerprinterTest.bitFingerprintTest.size();
        Assertions.assertEquals(64, tmpVariamycinFingerprintSize);
    }
    //
    /**
     * Tests whether the correct positions in the BitSet are set to true.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void bitSetTest() {
        BitSet tmpBitSetTest = new BitSet();
        tmpBitSetTest.set(3);
        tmpBitSetTest.set(5);
        tmpBitSetTest.set(9);
        tmpBitSetTest.set(14);
        tmpBitSetTest.set(16);
        tmpBitSetTest.set(17);
        tmpBitSetTest.set(18);
        tmpBitSetTest.set(26);
        tmpBitSetTest.set(27);
        BitSet tmpVariamycinBitSet = FragmentFingerprinterTest.bitFingerprintTest.asBitSet();
        Assertions.assertEquals(tmpBitSetTest, tmpVariamycinBitSet);
    }
    //
    /**
     * Tests whether all positively set positions are actually returned.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getBitSetTest() {
        int[] tmpArrayBitSetTest = {3,5,9,14,16,17,18,26,27};
        int[] tmpVariamcyinArrayBitSet = FragmentFingerprinterTest.bitFingerprintTest.getSetbits();
        Assertions.assertArrayEquals(tmpArrayBitSetTest, tmpVariamcyinArrayBitSet);
    }
    //
    /**
     * Tests the method getBitArray(List<String>)
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getBitArrayTestInputList() {
        int[] tmpTestArray = new int[28];
        tmpTestArray[3] = 1;
        tmpTestArray[5] = 1;
        tmpTestArray[9] = 1;
        tmpTestArray[14] = 1;
        tmpTestArray[16] = 1;
        tmpTestArray[17] = 1;
        tmpTestArray[18] = 1;
        tmpTestArray[26] = 1;
        tmpTestArray[27] = 1;
        int[] tmpVariamycinBitArray = FragmentFingerprinterTest.fragmentFingerprinter.getBitArray(FragmentFingerprinterTest.dataForGeneratingBitFingerprint);
        Assertions.assertArrayEquals(tmpTestArray, tmpVariamycinBitArray);
    }
    //
    /**
     * Tests the method getBitArray(Map<String,Integer>)
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getBitArrayTestInputMap() {
        int[] tmpTestArray = new int[28];
        tmpTestArray[3] = 1;
        tmpTestArray[5] = 1;
        tmpTestArray[9] = 1;
        tmpTestArray[14] = 1;
        tmpTestArray[16] = 1;
        tmpTestArray[17] = 1;
        tmpTestArray[18] = 1;
        tmpTestArray[26] = 1;
        tmpTestArray[27] = 1;
        int[] tmpVariamycinBitArray = FragmentFingerprinterTest.fragmentFingerprinter.getBitArray(FragmentFingerprinterTest.moleculeFragmentList.get(FragmentFingerprinterTest.moleculeFragmentList.size() - 1));
        Assertions.assertArrayEquals(tmpTestArray, tmpVariamycinBitArray);
    }
    //
    /**
     * Tests the size of the count fingerprint
     *
     * Test molecule: Variamycin
     *
     */
    @Test
    public void countFingerprintSizeTest() {
        long tmpCountFingerprintSizeTest = 28;
        long tmpVariamycinCountFingerprintSize = FragmentFingerprinterTest.countFingerprintTest.size();
        Assertions.assertEquals(tmpCountFingerprintSizeTest, tmpVariamycinCountFingerprintSize);
    }
    //
    /**
     * Tests the method numberOfPopulatedBins()
     *
     * Test molecule: Variamycin
     */
    @Test
    public void numberOfPopulatedBinsTest() {
        int tmpBinsTest = 28;
        int tmpVariamycinCountFingerprintNumberOfPopulatedBins = FragmentFingerprinterTest.countFingerprintTest.numOfPopulatedbins();
        Assertions.assertEquals(tmpBinsTest, tmpVariamycinCountFingerprintNumberOfPopulatedBins);
    }
    //
    /**
     * Tests the count value at position 17 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountTest() {
        int tmpCountTestForGivenIndex = 8;
        int tmpCountForGivenIndexInVariamycinFingerprint = FragmentFingerprinterTest.countFingerprintTest.getCount(17);
        Assertions.assertEquals(tmpCountTestForGivenIndex, tmpCountForGivenIndexInVariamycinFingerprint);
    }
    //
    /**
     * Tests the hash value at position 17 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getHashTest(){
        int tmpHashTestForGivenIndex = 17;
        int tmpHashForGivenIndexInVariamycin =  FragmentFingerprinterTest.countFingerprintTest.getHash(17);
        Assertions.assertEquals(tmpHashTestForGivenIndex, tmpHashForGivenIndexInVariamycin);
    }
    //
    /**
     * Tests whether the correct count value is supplied for the hash value 10.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountForHashTest() {
        int tmpCountForGivenHashValue = 0;
        int tmpCountForHashInVariamycinFingerprint = FragmentFingerprinterTest.countFingerprintTest.getCountForHash(10);
        Assertions.assertEquals(tmpCountForGivenHashValue, tmpCountForHashInVariamycinFingerprint);
    }
    //
    /**
     * Tests whether the fingerprint contains the given hash. The given hash is 30.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void hasHashTest() {
        boolean tmpHasHashForGivenIndex = false;
        boolean tmpHasHashForGivenIndexInVariamycinFingerprint = FragmentFingerprinterTest.countFingerprintTest.hasHash(30);
        Assertions.assertEquals(tmpHasHashForGivenIndex, tmpHasHashForGivenIndexInVariamycinFingerprint);
    }
    //
    /**
     * Tests the count method
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countTest() {
        int tmpCountForGivenSmilesString = 5;
        int tmpCountForGivenSmilesStringInVariamycinFingerprint = FragmentFingerprinterTest.countFingerprintTest.count("CCCCC");
        Assertions.assertEquals(tmpCountForGivenSmilesString, tmpCountForGivenSmilesStringInVariamycinFingerprint);
    }
    //
    /**
     * Tests the size of the fingerprint
     *
     * Test molecule: Variamycin
     */
    @Test
    public void fragmentFingerprintSizeTest() {
        int tmpFingerprintSizeTest = 28;
        int tmpVariamycinFingerprintSize = FragmentFingerprinterTest.fragmentFingerprinter.getSize();
        Assertions.assertEquals(tmpFingerprintSizeTest, tmpVariamycinFingerprintSize);
    }
    //
    /**
     * Tests the method getCountFingerprint, the input must be a list
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintSizeTestInputList() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles);
        long tmpSizeTest = 28;
        Assertions.assertEquals(tmpSizeTest, tmpCountFingerprintInputList.size());
    }
    //
    /**
     * Tests the method numberOfPopulatedBins()
     *
     * Test molecule: Variamycin
     */
    @Test
    public void numberOfPopulatedBinsTestInputList() {
        int tmpNumberOfPopulatedBinsVariamycin = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles).numOfPopulatedbins();
        int tmpBinsTest = 28;
        Assertions.assertEquals(tmpBinsTest, tmpNumberOfPopulatedBinsVariamycin);
    }
    //
    /**
     * Tests the count value at position 26 in the count fingerprint.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountTestInputList() {
        int tmpCountForGivenIndexTest = 5;
        int tmpCountForGivenIndexInVariamycinFingerprint = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles).getCount(26);
        Assertions.assertEquals(tmpCountForGivenIndexTest, tmpCountForGivenIndexInVariamycinFingerprint);
    }
    //
    /**
     * Tests the hash value at position 10 in the count fingerprint.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getHashTestInputList() {
        int tmpHashForGivenIndexTest = 10;
        int tmpHashForGivenIndexInVariamcyinFingerprint = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles).getHash(10);
        Assertions.assertEquals(tmpHashForGivenIndexTest,tmpHashForGivenIndexInVariamcyinFingerprint);
    }
    //
    /**
     * Tests whether the fingerprint contains the given hash. The given hash is 20.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void hasHashTestInputList() {
        boolean tmpHasHashForGivenIndex = true;
        boolean tmpHasHashForGivenIndexInVariamycinFingerprint = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles).hasHash(20);
        Assertions.assertEquals(tmpHasHashForGivenIndex, tmpHasHashForGivenIndexInVariamycinFingerprint);
    }
    //
    /**
     * Tests whether the correct count value is supplied for the hash value 0.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountForHashTestInputList() {
        int tmpCountForHashTest = 0;
        int tmpCountForHashInVariamycinFingerprint = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles).getCountForHash(0);
        Assertions.assertEquals(tmpCountForHashTest,tmpCountForHashInVariamycinFingerprint);
    }
    //
    /**
     * Tests the method getCountArray(List<String>).
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountArrayInputListTest() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[3] = 1;
        tmpTestCountArray[5] = 2;
        tmpTestCountArray[9] = 1;
        tmpTestCountArray[14] = 5;
        tmpTestCountArray[16] = 2;
        tmpTestCountArray[17] = 8;
        tmpTestCountArray[18] = 2;
        tmpTestCountArray[26] = 5;
        tmpTestCountArray[27] = 2;
        int[] tmpVariamycinCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.countListOfUniqueSmiles);
        Assertions.assertArrayEquals(tmpTestCountArray, tmpVariamycinCountArray);
    }
    //
    /**
     * Tests the method getCountArray(Map<String,Integer>)
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountArrayInputMapTest() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[3] = 1;
        tmpTestCountArray[5] = 2;
        tmpTestCountArray[9] = 1;
        tmpTestCountArray[14] = 5;
        tmpTestCountArray[16] = 2;
        tmpTestCountArray[17] = 8;
        tmpTestCountArray[18] = 2;
        tmpTestCountArray[26] = 5;
        tmpTestCountArray[27] = 2;
        int[] tmpVariamycinCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(FragmentFingerprinterTest.moleculeFragmentList.size() - 1));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpVariamycinCountArray);
    }
    //
    /**
     * Tests the method getBitDefinition()
     *
     */
    @Test
    public void getBitDefinitionTest() {
        String tmpBitDefinitionForGivenBitTest = "[H]Oc";
        String tmpBitDefinitionForGivenBitInVariamycinFingerprint = FragmentFingerprinterTest.fragmentFingerprinter.getBitDefinition(27);
        Assertions.assertEquals(tmpBitDefinitionForGivenBitTest, tmpBitDefinitionForGivenBitInVariamycinFingerprint);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Test count arrays of all molecules" defaultstate="collapsed">
    /**
     * Tests count array
     *
     * Test molecule: Valdiazen
     */
    @Test
    public void getValdiazenCountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[0] = 1;
        tmpTestCountArray[8] = 1;
        tmpTestCountArray[17] = 1;
        int[] tmpValdiazenCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(0));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpValdiazenCountArray);
    }
    //
    /**
     * Tests count array
     *
     * Test molecule: Napthomycin D
     */
    @Test
    public void getNapthomycinDCountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[15] = 1;
        tmpTestCountArray[25] = 2;
        tmpTestCountArray[17] = 3;
        tmpTestCountArray[27] = 1;
        tmpTestCountArray[5] = 3;
        tmpTestCountArray[7] = 1;
        tmpTestCountArray[21] = 1;
        tmpTestCountArray[12] = 1;
        tmpTestCountArray[6] = 1;
        tmpTestCountArray[26] = 1;
        int[] tmpNapthomycinDCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(1));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpNapthomycinDCountArray);
    }
    //
    /**
     * Tests count array
     *
     * Test molecule: Nona-2,6-dienal
     */
    @Test
    public void getNonaDienalCountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[21] = 1;
        tmpTestCountArray[12] = 2;
        tmpTestCountArray[19] = 1;
        int[] tmpNonaDienalCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(2));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpNonaDienalCountArray);
    }
    //
    /**
     * Tests count array
     *
     * Test molecule: Istanbulin A
     */
    @Test
    public void getIstanbulinACountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[16] = 1;
        tmpTestCountArray[13] = 1;
        tmpTestCountArray[5] = 1;
        tmpTestCountArray[24] = 1;
        int[] tmpIstanbulinACountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(3));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpIstanbulinACountArray);
    }
    //
    /**
     * Tests count array
     *
     * Test molecule: Estradiol
     */
    @Test
    public void getEstradiolCountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[17] = 1;
        tmpTestCountArray[27] = 1;
        tmpTestCountArray[4] = 1;
        int[] tmpEstradiolCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(4));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpEstradiolCountArray);
    }
    //
    /**
     * Tests count array
     *
     * Test molecule: Flower of Paradise
     */
    @Test
    public void getFlowerOfParadiseCountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[1] = 1;
        tmpTestCountArray[20] = 1;
        int[] tmpFlowerOfParadiseCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(5));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpFlowerOfParadiseCountArray);
    }
    //
    /**
     * Tests count array
     *
     * Test molecule: Curcumin
     */
    @Test
    public void getCurcuminCountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[25] = 2;
        tmpTestCountArray[18] = 2;
        tmpTestCountArray[27] = 2;
        tmpTestCountArray[5] = 3;
        tmpTestCountArray[20] = 2;
        int[] tmpCurcuminCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(6));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpCurcuminCountArray);
    }
    //
    /**
     * Tests count array
     *
     * Test molecule: Robinetidinol
     *
     */
    @Test
    public void getRobinetidinolCountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[18] = 2;
        tmpTestCountArray[17] = 2;
        tmpTestCountArray[27] = 8;
        tmpTestCountArray[10] = 1;
        int[] tmpRobinetidinolCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(7));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpRobinetidinolCountArray);
    }
    //
    /**
     * Tests count array
     *
     * Test molecule: Alkaloid
     */
    @Test
    public void getAlkaloidCountArray() {
        int[] tmpTestCountArray = new int[FragmentFingerprinterTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[23] = 1;
        tmpTestCountArray[2] = 1;
        tmpTestCountArray[11] = 1;
        tmpTestCountArray[5] = 1;
        tmpTestCountArray[3] = 1;
        tmpTestCountArray[22] = 1;
        int[] tmpAlkaloidCountArray = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(8));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpAlkaloidCountArray);
    }
    //</editor-fold>
    //
    //<editor-fold desc="Test bit and count fingerprint of naphthalene fingerprint" defaultstate="collapsed">
    //
    /**
     * Tests the size of the naphthalene fingerprint
     */
    @Test
    public void getNaphthaleneFingerprintSize() {
        int tmpNaphthaleneFingerprintSizeTest = 7;
        int tmpNaphthaleneFingerprintSize = FragmentFingerprinterTest.naphthaleneFingerprinter.getSize();
        Assertions.assertEquals(tmpNaphthaleneFingerprintSizeTest, tmpNaphthaleneFingerprintSize);
    }
    //
    /**
     * Tests the number of positive bits in the naphthalene fingerprint.
     */
    @Test
    public void getNumberOfPositiveBitsInNaphthaleneFingerprint() {
        int tmpNumberOfPositiveBitsTest = 2;
        int tmpNumberOfPositiveBits = FragmentFingerprinterTest.cNP0437667BitFP.cardinality();
        Assertions.assertEquals(tmpNumberOfPositiveBitsTest, tmpNumberOfPositiveBits);
    }
    //
    /**
     * Tests the bit fingerprint of the naphthalene derivate
     */
    @Test
    public void getNaphthaleneBitFingerprint() {
        int[] tmpBitFingerprintTest = new int[7];
        tmpBitFingerprintTest[0] = 0;
        tmpBitFingerprintTest[1] = 1;
        tmpBitFingerprintTest[2] = 0;
        tmpBitFingerprintTest[3] = 0;
        tmpBitFingerprintTest[4] = 1;
        tmpBitFingerprintTest[5] = 0;
        tmpBitFingerprintTest[6] = 0;
        int[] tmpBitFingerprint = FragmentFingerprinterTest.naphthaleneFingerprinter.getBitArray(FragmentFingerprinterTest.cNP0437667Fragments);
        Assertions.assertArrayEquals(tmpBitFingerprintTest, tmpBitFingerprint);
    }
    //</editor-fold>
}

