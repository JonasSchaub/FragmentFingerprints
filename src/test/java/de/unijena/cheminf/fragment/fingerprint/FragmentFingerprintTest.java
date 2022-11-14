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
        BufferedReader tmpFragmentSetReader = new BufferedReader(new FileReader("src/main/resources/de/unijena/cheminf/fragment/fingerprint/FragmentsList.txt"));
        String tmpLine;
        String tmpSeparatorComma = ",";
        while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
            String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
            FragmentFingerprintTest.fragmentList.add(tmpSmilesOfFragments[0]);
        }
        FragmentFingerprintTest.fragmentList.remove(0);
        // Read CSV file (itemization tab)
        BufferedReader tmpMoleculeFragmentsReader = new BufferedReader(new FileReader("src/main/resources/de/unijena/cheminf/fragment/fingerprint/MoleculeFragments.txt"));
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpList = new ArrayList<>();
        String tmpMoleculeLine;
        System.out.println("\n\tBit and count arrays of the given molecules:");
        while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
            String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
            tmpList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
        }
        List<String> tmpSeparateList;
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
            FragmentFingerprinter tmpFingerprintRepresentation = new FragmentFingerprinter(FragmentFingerprintTest.fragmentList);
            IBitFingerprint tmpBitFingerprint = tmpFingerprintRepresentation.getBitFingerprint(FragmentFingerprintTest.dataForGenerateBitFingerprint);
            ICountFingerprint tmpCountFingerprint = tmpFingerprintRepresentation.getCountFingerprint(tmpMoleculeFragmentsMap);
            System.out.println("\t\tBit array: "+ tmpSeparateList.get(0)+ ": "+ java.util.Arrays.toString(tmpFingerprintRepresentation.getBitArray(FragmentFingerprintTest.dataForGenerateBitFingerprint)));
            System.out.println("\t\tCount array: "+ tmpSeparateList.get(0)+ ": "+ java.util.Arrays.toString(tmpFingerprintRepresentation.getCountArray(tmpMoleculeFragmentsMap)));
            System.out.println("\t\t     ");

            // add all molecule maps
            FragmentFingerprintTest.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
        }
        // Objects necessary for the test are created (used only in @Test)
        FragmentFingerprintTest.fragmentFingerprinter = new FragmentFingerprinter(FragmentFingerprintTest.fragmentList);
        FragmentFingerprintTest.countFingerprintTest = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.moleculeFragmentList.get(FragmentFingerprintTest.moleculeFragmentList.size()-1));
        FragmentFingerprintTest.bitFingerprintTest = FragmentFingerprintTest.fragmentFingerprinter.getBitFingerprint(FragmentFingerprintTest.dataForGenerateBitFingerprint);
    }

    /**
     * Tests the getBitVectorTest() method
     * Test molecule: Variamycin
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void getBitVectorTest() throws Exception{
        // Test complete bit vector
        int[] tmpBitArrayTest = {0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
        Assert.assertArrayEquals(tmpBitArrayTest, FragmentFingerprintTest.fragmentFingerprinter.getBitArray(FragmentFingerprintTest.dataForGenerateBitFingerprint));
    }

    /**
     * Tests the getCountVectorTest() method, with a Map as input
     * Test molecule: Variamycin
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void getCountVectorTest() throws Exception{
        // Test complete count vector; Input is a map
        int[] tmpCountArrayTest = {0, 2, 0, 0, 0, 0, 2, 0, 0, 5, 8, 2, 0, 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 5, 0};
        Assert.assertArrayEquals(tmpCountArrayTest, FragmentFingerprintTest.fragmentFingerprinter.getCountArray(FragmentFingerprintTest.moleculeFragmentList.get(FragmentFingerprintTest.moleculeFragmentList.size()-1)));
    }
    /**
     * Tests the other methods used in the creation of the bit fingerprint
     * Test molecule: Variamycin
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void bitFingerprintTest() throws  Exception {
        //Test number of positive indices
        int tmpNumberPositiveBitsTest = 9;
        Assert.assertEquals(tmpNumberPositiveBitsTest, FragmentFingerprintTest.bitFingerprintTest.cardinality());
        // Test size
        long tmpSizeTest = 64;
        Assert.assertEquals(64,FragmentFingerprintTest.bitFingerprintTest.size());
        //Test asBitSet
        BitSet tmpBitSetTest = new BitSet();
        tmpBitSetTest.set(1);
        tmpBitSetTest.set(6);
        tmpBitSetTest.set(9);
        tmpBitSetTest.set(10);
        tmpBitSetTest.set(11);
        tmpBitSetTest.set(13);
        tmpBitSetTest.set(14);
        tmpBitSetTest.set(21);
        tmpBitSetTest.set(26);
        Assert.assertEquals(tmpBitSetTest,FragmentFingerprintTest.bitFingerprintTest.asBitSet());
        // Test getSetbits
        int[] tmpArrayBitSetTest = {1,6,9,10,11,13,14,21,26};
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
        // Test value at position 10
        Assert.assertEquals(5, FragmentFingerprintTest.countFingerprintTest.getCount(26));
        // Test the hash in index 10
        Assert.assertEquals(10,FragmentFingerprintTest.countFingerprintTest.getHash(10));
        // Test whether the fingerprint contains the given hash. index 0-27 true; index>27 false
        Assert.assertEquals(false, FragmentFingerprintTest.countFingerprintTest.hasHash(30));
        // Test the count value for the bin with index 10.
        Assert.assertEquals(0, FragmentFingerprintTest.countFingerprintTest.getCountForHash(0));
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
        Assert.assertEquals(FragmentFingerprintTest.fragmentList.size(),FragmentFingerprintTest.fragmentFingerprinter.getSize());
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
        tmpCountListOfUniqueSmiles.add("*C(*)=O");
        tmpCountListOfUniqueSmiles.add("*C(*)=O");
        tmpCountListOfUniqueSmiles.add("[H]Oc");
        tmpCountListOfUniqueSmiles.add("[H]Oc");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("*O*");
        tmpCountListOfUniqueSmiles.add("*O*");
        tmpCountListOfUniqueSmiles.add("c1cc(cc2ccc(cc12)CC(C)C)C");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("C");
        tmpCountListOfUniqueSmiles.add("C");
        tmpCountListOfUniqueSmiles.add("CCC");

        //Object
        ICountFingerprint tmpCountFingerprintInputList = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(tmpCountListOfUniqueSmiles);

        // Test size
        long tmpSizeTest = FragmentFingerprintTest.fragmentList.size();
        Assert.assertEquals(tmpSizeTest, tmpCountFingerprintInputList.size());
        // Test the number of bins that are populated.
        int tmpBinsTest = 28;
        Assert.assertEquals(tmpBinsTest, tmpCountFingerprintInputList.numOfPopulatedbins());
        // Test value at position 10
        Assert.assertEquals(5, tmpCountFingerprintInputList.getCount(26));
        // Test the hash in index 10
        Assert.assertEquals(10,tmpCountFingerprintInputList.getHash(10));
        // Test whether the fingerprint contains the given hash. index 0-27 true; index>27 false
        Assert.assertEquals(false, tmpCountFingerprintInputList.hasHash(30));
        // Test the count value for the bin with index 10.
        Assert.assertEquals(0, tmpCountFingerprintInputList.getCountForHash(0));
    }

    /**
     * Tests the method getCountArray, the input must be a list
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void countArrayTestInputList() throws Exception {
        // list to generate input for the  "getCountFingerprint(List<String> aListToCountMoleculeFragmentSmiles)" method
        ArrayList<String> tmpCountListOfUniqueSmiles = new ArrayList<>();
        tmpCountListOfUniqueSmiles.add("*C(*)=O");
        tmpCountListOfUniqueSmiles.add("*C(*)=O");
        tmpCountListOfUniqueSmiles.add("[H]Oc");
        tmpCountListOfUniqueSmiles.add("[H]Oc");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("*OCO*");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("[H]OC");
        tmpCountListOfUniqueSmiles.add("*O*");
        tmpCountListOfUniqueSmiles.add("*O*");
        tmpCountListOfUniqueSmiles.add("c1cc(cc2ccc(cc12)CC(C)C)C");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("CCCCC");
        tmpCountListOfUniqueSmiles.add("C");
        tmpCountListOfUniqueSmiles.add("C");
        tmpCountListOfUniqueSmiles.add("CCC");

        int[] tmpCountArrayTestIfInputIsList = {0, 2, 0, 0, 0, 0, 2, 0, 0, 5, 8, 2, 0, 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 5, 0};
        Assert.assertArrayEquals(tmpCountArrayTestIfInputIsList, FragmentFingerprintTest.fragmentFingerprinter.getCountArray(tmpCountListOfUniqueSmiles));
    }


}

