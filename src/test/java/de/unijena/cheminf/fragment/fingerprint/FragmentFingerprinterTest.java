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
    private static final int initialCapacityValue = 200;
    //</editor-fold>
    //
    //<editor-fold desc="private static class variables" defaultstate="collapsed">
    /**
     * List in which all molecule fragments are stored that are read in from the CSV file.
     */
    private static ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>(FragmentFingerprinterTest.initialCapacityValue); // TODO local variable if only test the last molecule
    /**
     * Is a list that contains all fragments, the fingerprint is then generated based on these fragments.
     */
    private static ArrayList<String> fragmentList = new ArrayList<>(FragmentFingerprinterTest.initialCapacityValue);
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
    public static void setUp() throws IOException {
        BufferedReader tmpFragmentSetReader;
        BufferedReader tmpMoleculeFragmentsReader;
        try {
            tmpFragmentSetReader = new BufferedReader(new FileReader("src/test/resources/de/unijena/cheminf/fragment/fingerprint/FragmentList.txt")); //src/test/resources/de/unijena/cheminf/fragment/fingerprint/FragmentList.txt
            tmpMoleculeFragmentsReader = new BufferedReader(new FileReader("src/test/resources/de/unijena/cheminf/fragment/fingerprint/MoleculeFragments.txt"));
        } catch (IOException anException) {
            throw new IOException("File is not readable!");
        }
        /* Read CSV file to obtain fragments used to create the fingerprint. A fragment in the form of
        unique SMILES is displayed at the first position of each line. Only these unique
        SMILES from the file are stored in the list.
         */
        String tmpLine;
        String tmpSeparatorComma = ",";
        try {
            while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
                String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
                FragmentFingerprinterTest.fragmentList.add(tmpSmilesOfFragments[0]);
            }
            FragmentFingerprinterTest.fragmentList.remove(0);
            tmpFragmentSetReader.close();
        } catch (IOException anException) {
            throw new IOException("invalid fragment file. At least one line is not readable.");
        }
        // Read CSV file
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpList = new ArrayList<>(FragmentFingerprinterTest.initialCapacityValue);
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
        FragmentFingerprinter tmpFingerprintRepresentation = new FragmentFingerprinter(FragmentFingerprinterTest.fragmentList);
        String tmpWorkingPath = (new File("src/test/resources/de/unijena/cheminf/fragment/fingerprint/").getAbsoluteFile().getAbsolutePath()) + File.separator;
        new File("src/test/resources/de/unijena/cheminf/fragment/fingerprint/" + "/Fingerprints").mkdirs();
        File tmpResultsLogFile = new File(tmpWorkingPath + "/Fingerprints/" + FragmentFingerprinterTest.FINGERPRINTS_FILE_NAME  + ".txt");
        FileWriter tmpResultsLogFileWriter = new FileWriter(tmpResultsLogFile, false);
        PrintWriter tmpResultPrintWriter = new PrintWriter(tmpResultsLogFileWriter);
        int[] tmpBitArray;
        int io = 0;
        for (int tmpCurrentLine = 1; tmpCurrentLine < tmpList.size(); tmpCurrentLine++) {
            tmpSeparateList = tmpList.get(tmpCurrentLine);
            List<String> ListWithoutNameAndMoleculeSmiles = tmpSeparateList.subList(2, tmpSeparateList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>();
            FragmentFingerprinterTest.dataForGenerationBitFingerprint = new ArrayList<>(FragmentFingerprinterTest.initialCapacityValue);
            for (int i = 0; i < ListWithoutNameAndMoleculeSmiles.size(); i++) {
                if (i % 2 == 0) {
                    tmpMoleculeFragmentsMap.put(ListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(ListWithoutNameAndMoleculeSmiles.get(i + 1)));
                    FragmentFingerprinterTest.dataForGenerationBitFingerprint.add(ListWithoutNameAndMoleculeSmiles.get(i));
                }
            }
            // Illustration the results of the bit arrays for the specified molecules
            try {
                IBitFingerprint tmpBitFingerprint = tmpFingerprintRepresentation.getBitFingerprint(FragmentFingerprinterTest.dataForGenerationBitFingerprint);
                ICountFingerprint tmpCountFingerprint = tmpFingerprintRepresentation.getCountFingerprint(tmpMoleculeFragmentsMap);
                tmpBitArray = tmpFingerprintRepresentation.getBitArray(FragmentFingerprinterTest.dataForGenerationBitFingerprint);
                int tmpLength = java.util.Arrays.toString(tmpBitArray).length();
                tmpCountFingerprint.setBehaveAsBitFingerprint(false);
                System.out.println("\t\tNumber of positive bits " + tmpSeparateList.get(0) + ": " + tmpBitFingerprint.cardinality());
                System.out.println("\t\tIndices of positive bits " + tmpSeparateList.get(0) + ": " + tmpBitFingerprint.asBitSet().toString());
                System.out.println("\t\tCount for the bin with the index 5 " + tmpSeparateList.get(0) + ": " + tmpCountFingerprint.getCount(5));
                tmpResultPrintWriter.println(java.util.Arrays.toString(tmpBitArray).substring(1, tmpLength - 1));
                tmpResultPrintWriter.flush();
            } catch( IllegalArgumentException  anException) {
                throw anException;
            } catch (NullPointerException anException) {
                System.out.println("Generation of one fingerprint did not work.");
            }
            // add all molecule maps
            FragmentFingerprinterTest.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
        }
        // Objects necessary for the test are created (used only in @Test)
        FragmentFingerprinterTest.fragmentFingerprinter = new FragmentFingerprinter(FragmentFingerprinterTest.fragmentList);
        FragmentFingerprinterTest.countFingerprintTest = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.moleculeFragmentList.get(FragmentFingerprinterTest.moleculeFragmentList.size() - 1));
        FragmentFingerprinterTest.bitFingerprintTest = FragmentFingerprinterTest.fragmentFingerprinter.getBitFingerprint(FragmentFingerprinterTest.dataForGenerationBitFingerprint);
        // Variamycin fragments
        FragmentFingerprinterTest.countListOfUniqueSmiles = new ArrayList<>(FragmentFingerprinterTest.initialCapacityValue);
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
        Assertions.assertEquals(tmpNumberPositiveBitsTest, FragmentFingerprinterTest.bitFingerprintTest.cardinality());
    }
    //
    /**
     * Tests the size of a bit fingerprint.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void bitFingerprintSizeTest() {
        Assertions.assertEquals(64, FragmentFingerprinterTest.bitFingerprintTest.size());
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
        Assertions.assertEquals(tmpBitSetTest, FragmentFingerprinterTest.bitFingerprintTest.asBitSet());
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
        Assertions.assertArrayEquals(tmpArrayBitSetTest, FragmentFingerprinterTest.bitFingerprintTest.getSetbits());
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
        Assertions.assertArrayEquals(tmpTestArray, FragmentFingerprinterTest.fragmentFingerprinter.getBitArray(FragmentFingerprinterTest.dataForGenerationBitFingerprint));
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
        Assertions.assertArrayEquals(tmpTestArray, FragmentFingerprinterTest.fragmentFingerprinter.getBitArray(FragmentFingerprinterTest.moleculeFragmentList.get(FragmentFingerprinterTest.moleculeFragmentList.size() - 1)));
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
        long tmpSizeTest = FragmentFingerprinterTest.fragmentList.size();
        Assertions.assertEquals(tmpSizeTest, FragmentFingerprinterTest.countFingerprintTest.size());
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
        Assertions.assertEquals(tmpBinsTest, FragmentFingerprinterTest.countFingerprintTest.numOfPopulatedbins());
    }
    //
    /**
     * Tests the count value at position 17 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountTest() {
        Assertions.assertEquals(8, FragmentFingerprinterTest.countFingerprintTest.getCount(17));
    }
    //
    /**
     * Tests the hash value at position 17 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getHashTest(){
        Assertions.assertEquals(17, FragmentFingerprinterTest.countFingerprintTest.getHash(17));
    }
    //
    /**
     * Tests whether the correct count value is supplied for the hash value 10.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountForHashTest() {
        Assertions.assertEquals(0, FragmentFingerprinterTest.countFingerprintTest.getCountForHash(10));
    }
    //
    /**
     * Tests whether the fingerprint contains the given hash. The given hash is 30.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void hasHashTest() {
        Assertions.assertEquals(false, FragmentFingerprinterTest.countFingerprintTest.hasHash(30));
    }
    //
    /**
     * Tests the size of the fingerprint
     *
     * Test molecule: Variamycin
     */
    @Test
    public void fragmentFingerprintSizeTest() {
        Assertions.assertEquals(28, FragmentFingerprinterTest.fragmentFingerprinter.getSize());
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
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles);
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
    public void getCountTestInputList() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles);
        Assertions.assertEquals(5, tmpCountFingerprintInputList.getCount(26));
    }
    //
    /**
     * Tests the hash value at position 10 in the count fingerprints.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getHashTestInputList() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles);
        Assertions.assertEquals(10,tmpCountFingerprintInputList.getHash(10));
    }
    //
    /**
     * Tests whether the fingerprint contains the given hash. The given hash is 20.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void hasHashTestInputList() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles);
        Assertions.assertEquals(true, tmpCountFingerprintInputList.hasHash(20));
    }
    //
    /**
     * Tests whether the correct count value is supplied for the hash value 0.
     *
     * Test molecule: Variamycin
     */
    @Test
    public void getCountForHashTestInputList() {
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprinterTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprinterTest.countListOfUniqueSmiles);
        Assertions.assertEquals(0, tmpCountFingerprintInputList.getCountForHash(0));
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
        Assertions.assertArrayEquals(tmpTestCountArray, FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.countListOfUniqueSmiles));
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
        FragmentFingerprinterTest.fragmentFingerprinter.getCountArray(FragmentFingerprinterTest.moleculeFragmentList.get(FragmentFingerprinterTest.moleculeFragmentList.size() - 1));
    }
    //
    /**
     * Tests the method getBitDefinition()
     *
     */
    @Test
    public void getBitDefinitionTest() {
        Assertions.assertEquals("*OC(*)=O", FragmentFingerprinterTest.fragmentFingerprinter.getBitDefinition(2));
    }
    //</editor-fold>
}

