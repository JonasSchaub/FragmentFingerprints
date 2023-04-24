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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
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
     * Initial capacity of the lists in which the data for generating the fingerprints are stored.
     */
    private static final int INITIAL_CAPACITY_VALUE = 10;
    /**
     *  Initial capacity of the lists in which stores //TODO
     */
    private static final int INITIAL_CAPACITY_VALUE_FOR_FRAGMENT_LIST = 28;
    /**
     * Value for determining the initial capacity of cards
     */
    private static final int INITIAL_CAPACITY_VALUE_FOR_MAP = Math.round((4/3) + 1);
    //</editor-fold>
    //
    //<editor-fold desc="private static class variables" defaultstate="collapsed">
    /**
     * List in which all molecule fragments are stored that are read in from the CSV file.
     */
    private static ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>(FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE);
    /**
     * Is a list that contains all fragments, the fingerprint is then generated based on these fragments.
     */
    private static ArrayList<String> fragmentList = new ArrayList<>(FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE_FOR_FRAGMENT_LIST);
    /**
     * fragmentFingerprinter
     */
    private static FragmentFingerprinter fragmentFingerprinter;
    /**
     * Bit fingerprint
     */
    private static IBitFingerprint bitFingerprintTest;
    /**
     * Count fingerprint
     */
    private static CountFingerprint countFingerprintTest;
    /**
     * List only with unique SMILES and without frequencies
     */
   private static ArrayList<String> dataForGeneratingBitFingerprint;
    /**
     * List in which a fragment occurs more than once
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
     * @BeforeAll ensures that the setUp method is only executed once.
     *
     * @throws IOException is thrown if the constructor is unable to open a text file for logging occurred exceptions.
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
        String tmpLine;
        String tmpSeparatorComma = ",";
        while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
            String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
            FragmentFingerprinterTest.fragmentList.add(tmpSmilesOfFragments[0]);
        }
        // removing header line value
        FragmentFingerprinterTest.fragmentList.remove(0);
        // Read CSV file
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpListOfMoleculesFragmentsAndFrequenciesList = new ArrayList<>(FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE);
        String tmpMoleculeLine;
        System.out.println("\n\tBit and count arrays of the given molecules:");
        while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
            String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
            tmpListOfMoleculesFragmentsAndFrequenciesList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
        }
        List<String> tmpMoleculeFragmentsAndFrequenciesList;
        // Instance is only used to print fingerprints on the console.
        FragmentFingerprinter tmpFragmentFingerprintRepresentation = new FragmentFingerprinter(FragmentFingerprinterTest.fragmentList);
        String tmpFingerprintOutputPath = (new File("").getAbsoluteFile().getAbsolutePath()) + File.separator;
        new File(tmpFingerprintOutputPath + "Fingerprints").mkdirs();
        File tmpFingerprintResultFile = new File(tmpFingerprintOutputPath + "/Fingerprints/" + FragmentFingerprinterTest.FINGERPRINTS_FILE_NAME  + ".txt");
        FileWriter tmpFingerprintResultsFileWriter = new FileWriter(tmpFingerprintResultFile, false);
        PrintWriter tmpFingerprintResultPrintWriter = new PrintWriter(tmpFingerprintResultsFileWriter);
        int[] tmpBitArray;
        for (int tmpCurrentLineIndex = 1; tmpCurrentLineIndex < tmpListOfMoleculesFragmentsAndFrequenciesList.size(); tmpCurrentLineIndex++) {
            tmpMoleculeFragmentsAndFrequenciesList = tmpListOfMoleculesFragmentsAndFrequenciesList.get(tmpCurrentLineIndex);
            List<String> tmpListWithoutNameAndMoleculeSmiles = tmpMoleculeFragmentsAndFrequenciesList.subList(2, tmpMoleculeFragmentsAndFrequenciesList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>(tmpListWithoutNameAndMoleculeSmiles.size()*INITIAL_CAPACITY_VALUE_FOR_MAP);
            FragmentFingerprinterTest.dataForGeneratingBitFingerprint = new ArrayList<>(FragmentFingerprinterTest.INITIAL_CAPACITY_VALUE);
            for (int i = 0; i < tmpListWithoutNameAndMoleculeSmiles.size(); i++) {
                if (i % 2 == 0) {
                    tmpMoleculeFragmentsMap.put(tmpListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(tmpListWithoutNameAndMoleculeSmiles.get(i + 1)));
                    FragmentFingerprinterTest.dataForGeneratingBitFingerprint.add(tmpListWithoutNameAndMoleculeSmiles.get(i));
                }
            }
            //Illustration of the results of the bit arrays for the specified molecules
            IBitFingerprint tmpBitFingerprint = tmpFragmentFingerprintRepresentation.getBitFingerprint(FragmentFingerprinterTest.dataForGeneratingBitFingerprint);
            ICountFingerprint tmpCountFingerprint = tmpFragmentFingerprintRepresentation.getCountFingerprint(tmpMoleculeFragmentsMap);
            tmpBitArray = tmpFragmentFingerprintRepresentation.getCountArray(FragmentFingerprinterTest.dataForGeneratingBitFingerprint);
            int tmpLength = java.util.Arrays.toString(tmpBitArray).length();
            tmpCountFingerprint.setBehaveAsBitFingerprint(false);
            System.out.println("\t\tNumber of positive bits " + tmpMoleculeFragmentsAndFrequenciesList.get(0) + ": " + tmpBitFingerprint.cardinality());
            System.out.println("\t\tIndices of positive bits " + tmpMoleculeFragmentsAndFrequenciesList.get(0) + ": " + tmpBitFingerprint.asBitSet().toString());
            System.out.println("\t\thas hash for the bin with the index 28 " + tmpMoleculeFragmentsAndFrequenciesList.get(0) + ": " + tmpCountFingerprint.hasHash(28));
           // System.out.println("\t\tHash for the bin with the index 28 " + tmpSeparateList.get(0) + ": " + tmpCountFingerprint.getHash(28));
            tmpFingerprintResultPrintWriter.println(java.util.Arrays.toString(tmpBitArray).substring(1, tmpLength - 1));
            tmpFingerprintResultPrintWriter.flush();
            // add all molecule maps
            FragmentFingerprinterTest.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
        }
        // Objects necessary for the test are created (used only in @Test)
        FragmentFingerprinterTest.fragmentFingerprinter = new FragmentFingerprinter(FragmentFingerprinterTest.fragmentList);
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
        System.out.println(bitFingerprintTest.cardinality()+"----cardinality tesst");
        System.out.println(fragmentFingerprinter.getVersionDescription() + "----version");

        System.out.println(java.util.Arrays.toString(fragmentFingerprinter.getCountArray(countListOfUniqueSmiles))+"----count ARray liste");
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
        long tmpFingerprintSizeOfVariamcyin =  FragmentFingerprinterTest.bitFingerprintTest.size();
        Assertions.assertEquals(64, tmpFingerprintSizeOfVariamcyin);
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
        BitSet tmpBitSetOfVariamycin = FragmentFingerprinterTest.bitFingerprintTest.asBitSet();
        Assertions.assertEquals(tmpBitSetTest,tmpBitSetOfVariamycin);
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
        int[] tmpVariamycinArray = FragmentFingerprinterTest.fragmentFingerprinter.getBitArray(FragmentFingerprinterTest.dataForGeneratingBitFingerprint);
        Assertions.assertArrayEquals(tmpTestArray, tmpVariamycinArray);
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
        long tmpCountFingerprintSizeOfVariamycin = FragmentFingerprinterTest.countFingerprintTest.size();
        Assertions.assertEquals(tmpCountFingerprintSizeTest, FragmentFingerprinterTest.countFingerprintTest.size());
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
        int tmpCountFingerprintNumberOfPopulatedBinsOfVariamycin = FragmentFingerprinterTest.countFingerprintTest.numOfPopulatedbins();
        Assertions.assertEquals(tmpBinsTest, tmpCountFingerprintNumberOfPopulatedBinsOfVariamycin);
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
        Assertions.assertEquals(false, tmpHasHashForGivenIndexInVariamycinFingerprint);
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
        int tmpFingerprintSizeOfVariamycin = FragmentFingerprinterTest.fragmentFingerprinter.getSize();
        Assertions.assertEquals(tmpFingerprintSizeTest, tmpFingerprintSizeOfVariamycin);
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
     * Tests the count value at position 26 in the count fingerprints.
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
     * Tests the hash value at position 10 in the count fingerprints.
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
        Assertions.assertEquals(true, tmpHasHashForGivenIndexInVariamycinFingerprint);
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
        Assertions.assertEquals(0,tmpCountForHashInVariamycinFingerprint);
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
        int[] tmpCountArrayOfVariamycin = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.countListOfUniqueSmiles);
        Assertions.assertArrayEquals(tmpTestCountArray, tmpCountArrayOfVariamycin);
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
        int[] tmpCountArrayOfVariamycin = FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(FragmentFingerprinterTest.moleculeFragmentList.size() - 1));
        Assertions.assertArrayEquals(tmpTestCountArray, tmpCountArrayOfVariamycin);
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
}

