/*
 * MIT License
 *
 * Copyright (c) 2022 Betuel Sevindik, Felix Baensch, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
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
package de.unijena.cheminf.fragment.fingerprint.main;

import de.unijena.cheminf.fragment.fingerprint.BitFragmentFingerprint;
import de.unijena.cheminf.fragment.fingerprint.CountFragmentFingerprint;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * Main class to start the application
 *
 * @author Betuel Sevindik
 */
public class Main {

    public static void main(String[] args) throws IOException {

        BitFragmentFingerprint tmpBitFragmentFingerprint = new BitFragmentFingerprint();
        CountFragmentFingerprint tmpCountFragmentFingerprint = new CountFragmentFingerprint();

        // read in a csv file created in MORTAR to get molecule fragments
        //Attention: the separator must be a semicolon
        BufferedReader tmpBufferedReader = new BufferedReader(new FileReader("C:\\Users\\betue\\Desktop\\Items_Ertl_algorithm_test.csv"));
        String tmpSeparator = ";";
        List<List<String>> tmpMoleculeFragments = new ArrayList<>();
        String tmpCsvLine;
        while ((tmpCsvLine = tmpBufferedReader.readLine()) != null) {
            String[] values = tmpCsvLine.split(tmpSeparator);
            tmpMoleculeFragments.add(Arrays.asList(values));
        }
        List<String> tmpMoleculeList = null;
        for(int tmpLine = 1; tmpLine< tmpMoleculeFragments.size(); tmpLine++) {
            tmpMoleculeList = tmpMoleculeFragments.get(tmpLine);
            //System.out.println(str + "---einezlene Listen");
            HashMap<String, String> tmpMoleculeMap = new HashMap<>();
            for(int tmpMoleculeListIndex = 0; tmpMoleculeListIndex < tmpMoleculeList.size(); tmpMoleculeListIndex++) {
                if (tmpMoleculeListIndex % 2 == 0) {
                    tmpMoleculeMap.put(tmpMoleculeList.get(tmpMoleculeListIndex), tmpMoleculeList.get(tmpMoleculeListIndex + 1));
                }
            }
           // System.out.println(allMaps +"------map");
            ArrayList tmpBitVector = tmpCountFragmentFingerprint.generateFragmentFingerprint((ArrayList<String>)getFragmentList(new File("C:\\Users\\betue\\Desktop\\Fragments_Ertl_algorithm.csv")), tmpMoleculeMap);
            System.out.println(tmpCountFragmentFingerprint.getIndicesPositiveBits(tmpBitVector)+ "----Positive indices");
            int[] arr = (int[]) tmpBitVector.get(0);
            System.out.println(java.util.Arrays.toString(arr) +"--"+ tmpMoleculeList.get(0));
            //System.out.println(java.util.Arrays.toString(g) +"--"+str.get(0));
           // int[] q = countFragmentFingerprint.generateFragmentFingerprint((ArrayList<String>) records, allMaps);
           // System.out.println(java.util.Arrays.toString(q)+ "---"+str.get(0));
        }

    }

    /**
     * Csv reader to get a list with fragments
     *
     * @param aFile
     * @return
     * @throws IOException
     */
    private static ArrayList<String> getFragmentList(File aFile) throws IOException {
        BufferedReader tmpBufferedReader = new BufferedReader(new FileReader(aFile));
        String tmpCsvLine = "";
        String DELIMITER = ",";
        ArrayList<String> tmpFragmentList = new ArrayList<>();
        while ((tmpCsvLine = tmpBufferedReader.readLine()) != null) {
            String[] values = tmpCsvLine.split(DELIMITER);
            tmpFragmentList.add(values[0]);
        }
        tmpFragmentList.remove(0);
        return tmpFragmentList;
    }
}
