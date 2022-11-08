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
import org.junit.Before;
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

    /**
     * List
     */
    ArrayList<HashMap<String, Integer>> moleculeSetList = new ArrayList<>();
    /**
     * Fragment list
     */
    ArrayList<String> fragmentList = new ArrayList<>();

    /**
     *  fragmentFingerprinter
     */
    FragmentFingerprinter fragmentFingerprinter;

    /**
     *  Bit fingerprint
     */
    IBitFingerprint bitFingerprintTest;

    /**
     *  Count fingerprint
     */
    ICountFingerprint countFingerprintTest;

    /**
     *  List only with unique SMILES and without frequencies
     */
    ArrayList<String> dataForGenerateBitFingerprint;

    /**
     * Empty Constructor
     */
    public FragmentFingerprintTest(){
    }

    /**
     * Start generating fingerprinting data (MORTAR)
     *
     * @throws IOException
     */
    @Before
    public void setUp() throws IOException {
        // Read CSV file ( fragmentation tab)
        BufferedReader tmpFragmentSetReader = new BufferedReader(new FileReader("src/main/resources/de/unijena/cheminf/fragment/fingerprint/FragmentsList.txt"));
        String tmpLine;
        String tmpSeparatorComma = ",";
        while ((tmpLine = tmpFragmentSetReader.readLine()) != null) {
            String[] tmpSmilesOfFragments = tmpLine.split(tmpSeparatorComma);
            fragmentList.add(tmpSmilesOfFragments[0]);
        }
        fragmentList.remove(0);
        // Read CSV file (itemization tab)
        BufferedReader tmpMoleculeFragmentsReader = new BufferedReader(new FileReader("src/main/resources/de/unijena/cheminf/fragment/fingerprint/MoleculeFragments.txt"));
        String tmpSeparatorSemicolon = ";";
        List<List<String>> tmpList = new ArrayList<>();
        String tmpMoleculeLine;
        while ((tmpMoleculeLine = tmpMoleculeFragmentsReader.readLine()) != null) {
            String[] tmpMoleculeFragmentsAndFrequencies = tmpMoleculeLine.split(tmpSeparatorSemicolon);
            tmpList.add(Arrays.asList(tmpMoleculeFragmentsAndFrequencies));
        }
        List<String> tmpSeparateList;
        for (int tmpCurrentLine = 1; tmpCurrentLine < tmpList.size(); tmpCurrentLine++) {
            tmpSeparateList = tmpList.get(tmpCurrentLine);
            List<String> ListWithoutNameAndMoleculeSmiles = tmpSeparateList.subList(2, tmpSeparateList.size());
            HashMap<String, Integer> tmpMoleculeFragmentsMap = new HashMap<>();
            this.dataForGenerateBitFingerprint = new ArrayList<>();
            for (int i = 0; i < ListWithoutNameAndMoleculeSmiles.size(); i++) {
                if (i % 2 == 0) {
                    tmpMoleculeFragmentsMap.put(ListWithoutNameAndMoleculeSmiles.get(i), Integer.valueOf(ListWithoutNameAndMoleculeSmiles.get(i + 1)));
                    this.dataForGenerateBitFingerprint.add(ListWithoutNameAndMoleculeSmiles.get(i));
                }
                /**
                 FragmentFingerprinter bitTest = new FragmentFingerprinter(records);
                 IBitFingerprint bitFingerprint = bitTest.getBitFingerprint(allMaps);
                 System.out.println(bitFingerprint.size()+ "----Size");
                 System.out.println(java.util.Arrays.toString(bitTest.getBitVector()) + "BitVector");
                 System.out.println(bitFingerprint.asBitSet().toString());
                 */
            }
            this.moleculeSetList.add(tmpMoleculeFragmentsMap);
        }
        // Objects necessary for the test are created
        this.fragmentFingerprinter = new FragmentFingerprinter(this.fragmentList);
        this.countFingerprintTest = this.fragmentFingerprinter.getCountFingerprint(this.moleculeSetList.get(this.moleculeSetList.size()-1));
        this.bitFingerprintTest = this.fragmentFingerprinter.getBitFingerprint(this.dataForGenerateBitFingerprint);
    }

    /**
     * Test the getBitVectorTest() method
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void getBitVectorTest() throws Exception{
        /**
            int it = 0;
            for (HashMap<String, Integer> tmpMap : moleculeSetList) {
                FragmentFingerprinter bitTest = new FragmentFingerprinter(this.fragmentsList);
                IBitFingerprint bitFingerprint = bitTest.getBitFingerprint(tmpMap);
                System.out.println(bitFingerprint.size()+ "----Size");
                System.out.println(java.util.Arrays.toString(bitTest.getBitVector()) + "BitVector");
                System.out.println(bitFingerprint.asBitSet().toString());
                it++;
            }
         */
        // Test complete bit vector
            int[] test = {0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
            Assert.assertArrayEquals(test, fragmentFingerprinter.getBitVector());
        }

    /**
     * Test the getCountVectorTest() method
     *
     * @throws Exception if anything goes wrong
     */
    @Test
        public void getCountVectorTest() throws Exception{
        // Test complete count vector
            int[] test = {0, 2, 0, 0, 0, 0, 2, 0, 0, 5, 8, 2, 0, 2, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 5, 0};
            Assert.assertArrayEquals(test, fragmentFingerprinter.getCountVector());
        }

    /**
     * Tests the other methods used in the creation of the bit fingerprint
     *
     * @throws Exception if anything goes wrong
     */
    @Test
        public void bitFingerprintTest() throws  Exception {
            //Test number of positive indices
            int tmpNumberPositiveBitsTest = 9;
            Assert.assertEquals(tmpNumberPositiveBitsTest, this.bitFingerprintTest.cardinality());
            // Test size
            long tmpSizeTest = 64;
            Assert.assertEquals(64,this.bitFingerprintTest.size());
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
            Assert.assertEquals(tmpBitSetTest,this.bitFingerprintTest.asBitSet());
            // Test getSetbits
            int[] tmpArrayBitSetTest = {1,6,9,10,11,13,14,21,26};
            Assert.assertArrayEquals(tmpArrayBitSetTest, this.bitFingerprintTest.getSetbits());
        }

    /**
     * Tests the other methods used in the creation of the count fingerprint
     *
     * @throws Exception if anything goes wrong
     */
    @Test
        public void countFingerprintTest() throws  Exception {
        // Test size
        long tmpSizeTest = this.fragmentList.size();
        Assert.assertEquals(tmpSizeTest, this.countFingerprintTest.size());
        // Test the number of bins that are populated.
        int tmpBinsTest = 28;
        Assert.assertEquals(tmpBinsTest, this.countFingerprintTest.numOfPopulatedbins());
        // Test value at position 10
        Assert.assertEquals(8, this.countFingerprintTest.getCount(10));
        // Test the hash in index 10
        Assert.assertEquals(10,this.countFingerprintTest.getHash(10));
        // Test whether the fingerprint contains the given hash. index 0-27 true; index>27 false
        Assert.assertEquals(true, this.countFingerprintTest.hasHash(27));
        // Test the count value for the bin with index 10.
        Assert.assertEquals(8, this.countFingerprintTest.getCount(10));
        }

    /**
     * Test the size of the fingerprint
     *
     * @throws Exception if anything goes wrong
     */
    @Test
        public void fragmentFingerprintSizeTest() throws Exception {
        // Test size of the fingerprint
        Assert.assertEquals(this.fragmentList.size(), this.fragmentFingerprinter.getSize());
        }
    }

