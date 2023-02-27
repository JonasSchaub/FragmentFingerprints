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
 */
public class FragmentFingerprintTest {
    //<editor-fold desc="private static final class variables" defaultstate="collapsed">
    /**
     * Name of file for writing fingerprints results.
     */
    private static final String FINGERPRINTS_FILE_NAME = "Fingerprints";
    /**
     * Initial capacity of the lists.
     */
    private static final int initialCapacityValue = 200;
    //</editor-fold>
    //
    //<editor-fold desc="private static class variables" defaultstate="collapsed">
    /**
     * Is a list in which all molecule fragments are stored that are read in from the CSV file.
     */
    private static ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>(FragmentFingerprintTest.initialCapacityValue); // TODO local variable if only test the last molecule
    /**
     * Is a list that contains all fragments, the fingerprint is then generated based on these fragments.
     */
    private static ArrayList<String> fragmentList = new ArrayList<>(FragmentFingerprintTest.initialCapacityValue);
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
    private static ICountFingerprint countFingerprintTest;
    /**
     * List only with unique SMILES and without frequencies
     */
   private static ArrayList<String> dataForGenerationBitFingerprint;
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
    public FragmentFingerprintTest() {
    }
    //</editor-fold>
    //
    //<editor-fold desc="BeforeAll method" defaultstate="collapsed">
    /**
     * Start generating fingerprinting data.
     * To create the fingerprints, 2 text files are used here. One of these text files contains the predefined
     * fragments and the other text file contains the fragments associated with the molecules.
     * Structure of the text files can be seen in resources folder
     *
     * @BeforeAll ensures that the setUp method is only executed once.
     *
     * @throws IOException is thrown if the constructor is unable to open a text file for logging occurred exceptions.
     */
    @BeforeAll
    public static void setUp() throws IOException {
        BufferedReader tmpFragmentSetReader;
        BufferedReader tmpMoleculeFragmentsReader;
        try {
            tmpFragmentSetReader = new BufferedReader(new FileReader("src/test/resources/de/unijena/cheminf/fragment/fingerprint/FragmentList.txt"));
            tmpMoleculeFragmentsReader = new BufferedReader(new FileReader("src/test/resources/de/unijena/cheminf/fragment/fingerprint/MoleculeFragments.txt"));
        } catch (IOException anException) {
            throw new IOException("File is not readable!");
        }
        // Read CSV file ( fragmentation tab) to  obtain fragments used to create the fingerprint.
        String tmpLine;
        String tmpSeparatorComma = ",";
        try {
            while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
                String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
                FragmentFingerprintTest.fragmentList.add(tmpSmilesOfFragments[0]);
            }
            FragmentFingerprintTest.fragmentList.remove(0);
            tmpFragmentSetReader.close();
        } catch (IOException anException) {
            throw new IOException("invalid fragment file. At least one line is not readable.");
        }
        // Read CSV file (itemization tab)
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpList = new ArrayList<>(FragmentFingerprintTest.initialCapacityValue);
        String tmpMoleculeLine;
        System.out.println("\n\tBit and count arrays of the given molecules:");
        try {
            while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
                String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
                tmpList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
            }
            tmpMoleculeFragmentsReader.close();
        } catch (IOException anException) {
            throw new IOException("invalid molecule file. At least one line is not readable");
        }
        List<String> tmpSeparateList;
        FragmentFingerprinter tmpFingerprintRepresentation = new FragmentFingerprinter(FragmentFingerprintTest.fragmentList);
        String tmpWorkingPath = (new File("src/test/resources/de/unijena/cheminf/fragment/fingerprint/").getAbsoluteFile().getAbsolutePath()) + File.separator;
        new File("src/test/resources/de/unijena/cheminf/fragment/fingerprint/" + "/Fingerprints").mkdirs();
        File tmpResultsLogFile = new File(tmpWorkingPath + "/Fingerprints/" + FragmentFingerprintTest.FINGERPRINTS_FILE_NAME  + ".txt");
        FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, false);
        PrintWriter tmpResultPrintWriter = new PrintWriter(tmpResultsLogFileWriter);
        int[] tmpBitArray;
        for (int tmpCurrentLine = 1; tmpCurrentLine < tmpList.size(); tmpCurrentLine++) {
            tmpSeparateList = tmpList.get(tmpCurrentLine);
            List<String> ListWithoutNameAndMoleculeSmiles = tmpSeparateList.subList(2, tmpSeparateList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>();
            FragmentFingerprintTest.dataForGenerationBitFingerprint = new ArrayList<>(FragmentFingerprintTest.initialCapacityValue);
            for (int i = 0; i < ListWithoutNameAndMoleculeSmiles.size(); i++) {
                if (i % 2 == 0) {
                    tmpMoleculeFragmentsMap.put(ListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(ListWithoutNameAndMoleculeSmiles.get(i + 1)));
                    FragmentFingerprintTest.dataForGenerationBitFingerprint.add(ListWithoutNameAndMoleculeSmiles.get(i));
                }
            }
            // Illustration the results of the bit arrays for the specified molecules
            try {
                IBitFingerprint tmpBitFingerprint = tmpFingerprintRepresentation.getBitFingerprint(FragmentFingerprintTest.dataForGenerationBitFingerprint);
                ICountFingerprint tmpCountFingerprint = tmpFingerprintRepresentation.getCountFingerprint(tmpMoleculeFragmentsMap);
                tmpBitArray = tmpFingerprintRepresentation.getBitArray();
                int tmpLength = java.util.Arrays.toString(tmpBitArray).length();
                tmpCountFingerprint.setBehaveAsBitFingerprint(false);
                System.out.println("\t\tNumber of positive bits " + tmpSeparateList.get(0) + ": " + tmpBitFingerprint.cardinality());
                System.out.println("\t\tIndices of positive bits " + tmpSeparateList.get(0) + ": " + tmpBitFingerprint.asBitSet().toString());
                System.out.println("\t\tCount for the bin with the index 5 " + tmpSeparateList.get(0) + ": " + tmpCountFingerprint.getCount(5));
                System.out.println("\t\tget hash" + tmpSeparateList.get(0) + ": " + tmpCountFingerprint.getHash(5));
                tmpResultPrintWriter.println(java.util.Arrays.toString(tmpBitArray).substring(1, tmpLength - 1));
                tmpResultPrintWriter.flush();
            } catch( IllegalArgumentException  anException) {
                throw anException;
            } catch (NullPointerException anException) {
                System.out.println("Generation of one fingerprint did not work.");
            }
            // add all molecule maps
            FragmentFingerprintTest.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
        }
        // Objects necessary for the test are created (used only in @Test)
        FragmentFingerprintTest.fragmentFingerprinter = new FragmentFingerprinter(FragmentFingerprintTest.fragmentList);
        FragmentFingerprintTest.countFingerprintTest = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.moleculeFragmentList.get(FragmentFingerprintTest.moleculeFragmentList.size() - 1));
        FragmentFingerprintTest.bitFingerprintTest = FragmentFingerprintTest.fragmentFingerprinter.getBitFingerprint(FragmentFingerprintTest.dataForGenerationBitFingerprint);
        // Variamycin fragments
        FragmentFingerprintTest.countListOfUniqueSmiles = new ArrayList<>(FragmentFingerprintTest.initialCapacityValue);
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*C(*)=O");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*C(*)=O");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]Oc");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]Oc");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*OCO*");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*O*");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("c1cc(cc2ccc(cc12)CC(C)C)C");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]OC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("CCCCC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("C");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("*O*");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("C");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("CCC");
        FragmentFingerprintTest.countListOfUniqueSmiles.add("[H]OC");
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
    public void bitFingerprintNumberOfPositiveIndicesTest() {
        int tmpNumberPositiveBitsTest = 9;
        Assertions.assertEquals(tmpNumberPositiveBitsTest, FragmentFingerprintTest.bitFingerprintTest.cardinality());
    }
    //
    /**
     * Tests the size of a bit fingerprint.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void bitFingerprintSizeTest() {
        Assertions.assertEquals(64,FragmentFingerprintTest.bitFingerprintTest.size());
    }
    //
    /**
     * Tests whether the correct positions in the BitSet are set to true.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void bitFingerprintBitSetTest() {
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
        Assertions.assertEquals(tmpBitSetTest,FragmentFingerprintTest.bitFingerprintTest.asBitSet());
    }
    //
    /**
     * Tests whether all positively set positions are actually returned.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void bitFingerprintGetBitSetTest() {
        int[] tmpArrayBitSetTest = {3,5,9,14,16,17,18,26,27};
        Assertions.assertArrayEquals(tmpArrayBitSetTest, FragmentFingerprintTest.bitFingerprintTest.getSetbits());
    }
    //
    /**
     * Tests the generated bit array for correctness.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void bitFingerprintBitArrayTest() {
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
        Assertions.assertArrayEquals(tmpTestArray, FragmentFingerprintTest.fragmentFingerprinter.getBitArray());
    }
    //
    /**
     * Tests the other methods used in the creation of the count fingerprint, if the input is a map
     *
     * Test molecule: Variamycin
     *
     */
    @Test
    public void countFingerprintSizeTest() {
        long tmpSizeTest = FragmentFingerprintTest.fragmentList.size();
        Assertions.assertEquals(tmpSizeTest, FragmentFingerprintTest.countFingerprintTest.size());
    }
    //
    /**
     * Tests the method numberOfPopulatedBins()
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintNumberOfPopulatedBinsTest() {
        int tmpBinsTest = 28;
        Assertions.assertEquals(tmpBinsTest, FragmentFingerprintTest.countFingerprintTest.numOfPopulatedbins());
    }
    //
    /**
     * Tests the count value at position 17 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintCountValueTest() {
        Assertions.assertEquals(8, FragmentFingerprintTest.countFingerprintTest.getCount(17));
    }
    //
    /**
     * Tests the hash value at position 17 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintHashValueTest(){
        Assertions.assertEquals(17,FragmentFingerprintTest.countFingerprintTest.getHash(17));
    }
    //
    /**
     * Tests whether the correct count value is supplied for the hash value 10.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintGetCountForHashTest() {
        Assertions.assertEquals(0, FragmentFingerprintTest.countFingerprintTest.getCountForHash(10));
    }
    //
    /**
     * Tests whether the fingerprint contains the given hash. The given hash is 30.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintHasHashValueTest() {
        Assertions.assertEquals(false, FragmentFingerprintTest.countFingerprintTest.hasHash(30));
    }
    //
    /**
     * Tests the size of the fingerprint
     *
     * Test molecule: Variamycin
     */
    @Test
    public void fragmentFingerprintSizeTest() {
        Assertions.assertEquals(28,FragmentFingerprintTest.fragmentFingerprinter.getSize());
    }
    //
    /**
     * Tests the method getCountFingerprint, the input must be a list
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintInputListSizeTest() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.countListOfUniqueSmiles);
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
    public void countFingerprintInputListNumberOfPopulatedBins() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.countListOfUniqueSmiles);
        int tmpBinsTest = 28;
        Assertions.assertEquals(tmpBinsTest, tmpCountFingerprintInputList.numOfPopulatedbins());
    }
    //
    /**
     * Tests the count value at position 26 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintInputListCountValueTest() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.countListOfUniqueSmiles);
        Assertions.assertEquals(5, tmpCountFingerprintInputList.getCount(26));
    }
    //
    /**
     * Tests the hash value at position 10 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintInputListGetHashTest() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.countListOfUniqueSmiles);
        Assertions.assertEquals(10,tmpCountFingerprintInputList.getHash(10));
    }
    //
    /**
     * Tests whether the fingerprint contains the given hash. The given hash is 20.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintInputListHasHashTest() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.countListOfUniqueSmiles);
        Assertions.assertEquals(true, tmpCountFingerprintInputList.hasHash(20));
    }
    //
    /**
     * Tests whether the correct count value is supplied for the hash value 0.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void countFingerprintInputListGetCountForHashTest() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.countListOfUniqueSmiles);
        Assertions.assertEquals(0, tmpCountFingerprintInputList.getCountForHash(0));
    }
    //
    /**
     * Tests the method getBitArray().
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getBitArrayTest() {
        //Test bit array
        int[] tmpTestBitArray = new int[FragmentFingerprintTest.fragmentFingerprinter.getSize()];
        tmpTestBitArray[3] = 1;
        tmpTestBitArray[5] = 1;
        tmpTestBitArray[9] = 1;
        tmpTestBitArray[14] = 1;
        tmpTestBitArray[16] = 1;
        tmpTestBitArray[17] = 1;
        tmpTestBitArray[18] = 1;
        tmpTestBitArray[26] = 1;
        tmpTestBitArray[27] = 1;
        Assertions.assertArrayEquals(tmpTestBitArray,FragmentFingerprintTest.fragmentFingerprinter.getBitArray());
    }
    //
    /**
     * Tests the method getCountArray().
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountArrayTest() {
        int[] tmpTestCountArray = new int[FragmentFingerprintTest.fragmentFingerprinter.getSize()];
        tmpTestCountArray[3] = 1;
        tmpTestCountArray[5] = 2;
        tmpTestCountArray[9] = 1;
        tmpTestCountArray[14] = 5;
        tmpTestCountArray[16] = 2;
        tmpTestCountArray[17] = 8;
        tmpTestCountArray[18] = 2;
        tmpTestCountArray[26] = 5;
        tmpTestCountArray[27] = 2;
        Assertions.assertArrayEquals(tmpTestCountArray,FragmentFingerprintTest.fragmentFingerprinter.getCountArray());
    }
    //</editor-fold>
}

