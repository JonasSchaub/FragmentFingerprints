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

import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Array;
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
        ArrayList<String> liste = new ArrayList<>();
        liste.add("CCC(C)C");
        liste.add("*OC(*)=O");
        liste.add("C");
        liste.add("*n(*)*");
        liste.add("[H]Oc");
        liste.add("[H]OC");
        liste.add("CCCC");
        liste.add("C=C");
        liste.add("*O*");
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
            FragmentFingerprinter tmpFingerprintRepresentation = new FragmentFingerprinter(liste); //FragmentFingerprintTest.fragmentList
           // IBitFingerprint tmpBitFingerprint = tmpFingerprintRepresentation.getBitFingerprint(FragmentFingerprintTest.dataForGenerateBitFingerprint);
            ICountFingerprint tmpCountFingerprint = tmpFingerprintRepresentation.getCountFingerprint(tmpMoleculeFragmentsMap);
            System.out.println(tmpCountFingerprint.getCount(8));
            System.out.println(tmpCountFingerprint.getHash(3));
            System.out.println(tmpCountFingerprint.hasHash(8));
            System.out.println(tmpCountFingerprint.size());
            // add all molecule maps
            FragmentFingerprintTest.moleculeFragmentList.add(tmpMoleculeFragmentsMap);
        }
        // Objects necessary for the test are created (used only in @Test)
        /*
        FragmentFingerprintTest.fragmentFingerprinter = new FragmentFingerprinter(FragmentFingerprintTest.fragmentList); //FragmentFingerprintTest.fragmentList
        FragmentFingerprintTest.countFingerprintTest = FragmentFingerprintTest.fragmentFingerprinter.getCountFingerprint(FragmentFingerprintTest.moleculeFragmentList.get(FragmentFingerprintTest.moleculeFragmentList.size()-1));
        FragmentFingerprintTest.bitFingerprintTest = FragmentFingerprintTest.fragmentFingerprinter.getBitFingerprint(FragmentFingerprintTest.dataForGenerateBitFingerprint);

         */
    }
    /**
     * Tests the getBitVectorTest() method
     * Test molecule: Variamycin
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void getBitVectorTest() throws Exception{ //TODO rename
        // Test complete bit vector
        /*
        int[] tmpBitArrayTest = {0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};
         assertArrayEquals(tmpBitArrayTest, FragmentFingerprintTest.fragmentFingerprinter.getBitArray());
         */

    }
}

