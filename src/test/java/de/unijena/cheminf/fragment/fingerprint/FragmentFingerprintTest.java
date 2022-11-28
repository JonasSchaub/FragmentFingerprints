/*
 * MIT License
 *
 * Copyright (c) 2022
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

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;

/**
 *  Class to test the correct working of  FragmentFingerprinter
 */
public class FragmentFingerprintTest  {
    //<editor-fold desc="private static class variables" defaultstate="collapsed">
    /**
     * Is a list in which all molecule fragments are stored that are read in from the CSV file.
     */
    static ArrayList<HashMap<String, Integer>> moleculeFragmentList = new ArrayList<>(); // TODO local variable if only test the last molecule
    /**
     * Is a list that contains all fragments, the fingerprint is then generated based on these fragments.
     */
    static ArrayList<String> fragmentList = new ArrayList<>();

    /**
     *  fragmentFingerprinter
     */
    static FragmentFingerprinter fragmentFingerprinter;

    /**
     *  Bit fingerprint
     */
    static IBitFingerprint bitFingerprintTest;

    /**
     *  Count fingerprint
     */
    static ICountFingerprint countFingerprintTest;

    /**
     *  List only with unique SMILES and without frequencies
     */
   static ArrayList<String> dataForGenerateBitFingerprint;
    //</editor-fold>
    //
    /**
     * Empty Constructor
     */
    public FragmentFingerprintTest(){
    }

    /**
     * Start generating fingerprinting data (MORTAR); @BeforeClass ensures that the setUp method is only executed once.
     *
     * @throws IOException
     */
    @BeforeClass
    public static void setUp() throws IOException {
        // Read CSV file ( fragmentation tab) to  obtain fragments used to create the fingerprint.
        BufferedReader tmpFragmentSetReader = new BufferedReader(new FileReader("src/test/resources/de/unijena/cheminf/fragment/fingerprint/FragmentList.txt"));
        String tmpLine;
        String tmpSeparatorComma = ",";
        while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
            String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
            FragmentFingerprintTest.fragmentList.add(tmpSmilesOfFragments[0]);
        }
        FragmentFingerprintTest.fragmentList.remove(0);
        // Read CSV file (itemization tab)
        BufferedReader tmpMoleculeFragmentsReader = new BufferedReader(new FileReader("src/test/resources/de/unijena/cheminf/fragment/fingerprint/MoleculeFragments.txt"));
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpList = new ArrayList<>();
        String tmpMoleculeLine;
        System.out.println("\n\tBit and count arrays of the given molecules:");
        while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
            String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
            tmpList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
        }
        List<String> tmpSeparateList;
        FragmentFingerprinter tmpFingerprintRepresentation = new FragmentFingerprinter(FragmentFingerprintTest.fragmentList);
        for (int tmpCurrentLine = 1; tmpCurrentLine < tmpList.size(); tmpCurrentLine++) {
            tmpSeparateList = tmpList.get(tmpCurrentLine);
            List<String> ListWithoutNameAndMoleculeSmiles = tmpSeparateList.subList(2, tmpSeparateList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>();
            FragmentFingerprintTest.dataForGenerateBitFingerprint = new ArrayList<>();
            for (int i = 0; i < ListWithoutNameAndMoleculeSmiles.size(); i++) {
                if (i % 2 == 0) {
                    tmpMoleculeFragmentsMap.put(ListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(ListWithoutNameAndMoleculeSmiles.get(i + 1)));
                    FragmentFingerprintTest.dataForGenerateBitFingerprint.add(ListWithoutNameAndMoleculeSmiles.get(i));
                }
            }
            // Illustration the results of the bit arrays for the specified molecules
            IBitFingerprint tmpBitFingerprint = tmpFingerprintRepresentation.getBitFingerprint(FragmentFingerprintTest.dataForGenerateBitFingerprint);
            ICountFingerprint tmpCountFingerprint = tmpFingerprintRepresentation.getCountFingerprint(tmpMoleculeFragmentsMap);
            System.out.println("\t\tNumber of positive bits " + tmpSeparateList.get(0) + ": " + tmpBitFingerprint.cardinality());
            System.out.println("\t\tIndices of positive bits "+ tmpSeparateList.get(0) + ": "  + tmpBitFingerprint.asBitSet().toString());
            System.out.println("\t\tCount for the bin with the index 5 " + tmpSeparateList.get(0) + ": " + tmpCountFingerprint.getCountForHash(5));
            // add all molecule maps
            FragmentFingerprintTest.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
        }
        // Objects necessary for the test are created (used only in @Test)
        FragmentFingerprintTest.fragmentFingerprinter = new FragmentFingerprinter(FragmentFingerprintTest.fragmentList); //FragmentFingerprintTest.fragmentList
        FragmentFingerprintTest.countFingerprintTest = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.moleculeFragmentList.get(FragmentFingerprintTest.moleculeFragmentList.size()-1));
        FragmentFingerprintTest.bitFingerprintTest = FragmentFingerprintTest.fragmentFingerprinter.getBitFingerprint(FragmentFingerprintTest.dataForGenerateBitFingerprint);
    }
    /**
     * Tests the getBitFingerprint() method
     * Test molecule: Variamycin
     *
     * @throws Exception if anything goes wrong
     */
    @Test
   public void bitFingerprintTest() {
        //Test number of positive indices
        int tmpNumberPositiveBitsTest = 9;
        Assert.assertEquals(tmpNumberPositiveBitsTest, FragmentFingerprintTest.bitFingerprintTest.cardinality());
        // Test size
        Assert.assertEquals(64,FragmentFingerprintTest.bitFingerprintTest.size());
        //Test asBitSet
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
        Assert.assertEquals(tmpBitSetTest,FragmentFingerprintTest.bitFingerprintTest.asBitSet());
        // Test getSetbits
        int[] tmpArrayBitSetTest = {3,5,9,14,16,17,18,26,27};
        Assert.assertArrayEquals(tmpArrayBitSetTest, FragmentFingerprintTest.bitFingerprintTest.getSetbits());
    }
    /**
     * Tests the other methods used in the creation of the count fingerprint, if the input is a map
     * Test molecule: Variamycin
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void countFingerprintTest() throws  Exception {
        // Test size
        long tmpSizeTest = FragmentFingerprintTest.fragmentList.size();
        Assert.assertEquals(tmpSizeTest, FragmentFingerprintTest.countFingerprintTest.size());
        // Test the number of bins that are populated.
        int tmpBinsTest = 28;
        Assert.assertEquals(tmpBinsTest, FragmentFingerprintTest.countFingerprintTest.numOfPopulatedbins());
        // Test value in bin 17
        Assert.assertEquals(8, FragmentFingerprintTest.countFingerprintTest.getCount(17));
        // Test the hash in index 17
        Assert.assertEquals(17,FragmentFingerprintTest.countFingerprintTest.getHash(17));
        // Test whether the fingerprint contains the given hash. index 0-27 true; index>27 false
        Assert.assertEquals(false, FragmentFingerprintTest.countFingerprintTest.hasHash(30));
        // Test the count value for the bin with index 10.
        Assert.assertEquals(0, FragmentFingerprintTest.countFingerprintTest.getCountForHash(10));
    }
    /**
     * Tests the size of the fingerprint
     * Test molecule: Variamycin
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void fragmentFingerprintSizeTest() throws Exception {
        // Test size of the fingerprint
        Assert.assertEquals(28,FragmentFingerprintTest.fragmentFingerprinter.getSize());
    }
    /**
     * Tests the method getCountFingerprint, the input must be a list
     *
     * Test molecule: Variamycin
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void countFingerprintTestInputList() throws Exception {
        // list to generate input for the  "getCountFingerprint(List<String> aListToCountMoleculeFragmentSmiles)" method
        ArrayList<String> tmpCountListOfUniqueSmiles = new ArrayList<>();
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("*C(*)=O");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("*C(*)=O");
        tmpCountListOfUniqueSmiles.add("[H]Oc");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("[H]Oc");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("*O*");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("c1cc(cc2ccc(cc12)CC(C)C)C");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("C");
        tmpCountListOfUniqueSmiles.add("*O*");
        tmpCountListOfUniqueSmiles.add("C");
        tmpCountListOfUniqueSmiles.add("CCC");
        tmpCountListOfUniqueSmiles.add("[H]OC");

        //Object
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(tmpCountListOfUniqueSmiles);

        // Test size
        long tmpSizeTest = FragmentFingerprintTest.fragmentList.size();
        Assert.assertEquals(tmpSizeTest, tmpCountFingerprintInputList.size());
        // Test the number of bins that are populated.
        int tmpBinsTest = 28;
        Assert.assertEquals(tmpBinsTest, tmpCountFingerprintInputList.numOfPopulatedbins());
        // Test value at position 26
        Assert.assertEquals(5, tmpCountFingerprintInputList.getCount(26));
        // Test the hash in index 10
        Assert.assertEquals(10,tmpCountFingerprintInputList.getHash(10));
        // Test whether the fingerprint contains the given hash. index 0-27 true; index>27 false
        Assert.assertEquals(false, tmpCountFingerprintInputList.hasHash(30));
        // Test the count value for the bin with index 10.
        Assert.assertEquals(0, tmpCountFingerprintInputList.getCountForHash(0));
    }



}

