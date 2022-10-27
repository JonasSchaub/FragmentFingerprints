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
package de.unijena.cheminf.fragment.fingerprint;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/**
 * Class to generate a bit fragment fingerprint
 *
 * @author Betuel Sevindik
 */
public class BitFragmentFingerprint implements IFingerprint {

    //<editor-fold desc="Constructor">
    /**
     * Empty constructor
     */
    public BitFragmentFingerprint() {
    }
    //</editor-fold>

    //<editor-fold desc="Public methods">
    /**
     * Method to generate a bit fragment fingerprint
     *
     * @param aFragmentList  a fragment list  containing all fragments
     * @param aMoleculeFragments  a molecule HashMap which contain as key the molecule fragments and as
     * value the frequency of the according fragment
     * @return an ArrayList; index 0 stores the fingerprint in the
     * form of an array and index 1 stores the positive indices in a list.
     */
    @Override
    public ArrayList generateFragmentFingerprint(ArrayList<String> aFragmentList, HashMap<String, String> aMoleculeFragments) {
        int tmpVectorSize = aFragmentList.size();
        int[] tmpBitVector = new int[tmpVectorSize];
        ArrayList<Integer> tmpPositiveIndices = new ArrayList<>();
        ArrayList<Object> tmpFingerprintResult = new ArrayList<>();
        for (int tmpDefaultBit = 0; tmpDefaultBit < tmpVectorSize; tmpDefaultBit++) {
            tmpBitVector[tmpDefaultBit] = 0;
        }
        HashMap<String, Integer> tmpFragmentHashMap = new HashMap<>();
        // sort aFragmentList alphabetically
        Collections.sort(aFragmentList, String.CASE_INSENSITIVE_ORDER);
        int tmpValuePosition = 0;
        for (String tmpKey : aFragmentList) {
            tmpFragmentHashMap.put(tmpKey, tmpValuePosition);
            tmpValuePosition++;
        }
        for (String tmpUniqueSmiles : aMoleculeFragments.keySet()) {
            if (tmpFragmentHashMap.containsKey(tmpUniqueSmiles) == true) {
                int tmpPosition = tmpFragmentHashMap.get(tmpUniqueSmiles);
                tmpPositiveIndices.add(tmpPosition);
                tmpBitVector[tmpPosition] = 1;
            }
        }
        tmpFingerprintResult.add(tmpBitVector);
        tmpFingerprintResult.add(tmpPositiveIndices);
        return tmpFingerprintResult;
    }
    //
    /**
     *  TODO implementation
     * @return
     */
    @Override
    public int getFragmentFingerprintSize(int[] aVector) {
        return aVector.length;
    }
    //
    /**
     * 
     * @return
     */
    public int getNumberPositiveBits(ArrayList<Object> aList) {
        Object tmpPositiveBits = aList.get(1);
        ArrayList tmpNumber  = (ArrayList) tmpPositiveBits;
        return tmpNumber.size();
    }
    //
    /**
     * Returns only the positive indices in the fingerprint
     *
     * @param aList
     * @return
     */
    @Override
    public ArrayList getIndicesPositiveBits(ArrayList<Object> aList) {
        Object tmpPositiveBits = aList.get(1);
        ArrayList<Integer> tmpBitsList = (ArrayList) tmpPositiveBits;
        return tmpBitsList;
    }
    //
    /**
     * TODO implementation
     *
     * @return
     */
    @Override
    public double calculateTanimotoSimilarity(ArrayList<Object> aFirstMoleculeFingerprint, ArrayList<Object> aSecondMoleculeFingerprint) {
        ArrayList<Integer> tmpFirstFingerprint = (ArrayList) aFirstMoleculeFingerprint;
        ArrayList<Integer> tmpSecondFingerprint = (ArrayList) aSecondMoleculeFingerprint;
        System.out.println(tmpFirstFingerprint +"---liste1");
        System.out.println(tmpSecondFingerprint +"---liste2");
        HashMap<Integer, Integer> tmpMapForFirstFingerprint = new HashMap<>();
        HashMap<Integer, Integer> tmpMapForSecondFingerprint = new HashMap<>();
        int tmpFirstMapValues = 0;
        int tmpSecondMapValues = 0;
        for (int tmpFirstMapKey : tmpFirstFingerprint) {
            tmpMapForFirstFingerprint.put(tmpFirstMapKey, tmpFirstMapValues);
            tmpFirstMapValues++;
        }
        System.out.println(tmpMapForFirstFingerprint + "-----map1");
        for (int tmpSecondMapKey : tmpSecondFingerprint) {
            tmpMapForSecondFingerprint.put(tmpSecondMapKey, tmpSecondMapValues);
            tmpSecondMapValues++;
        }
        System.out.println(tmpMapForSecondFingerprint +"----map2");
        int tmpIntersection = 0;
            for (int tmpEqualsKey : tmpMapForSecondFingerprint.keySet()) {
                if (tmpMapForFirstFingerprint.containsKey(tmpEqualsKey)) {
                    tmpIntersection++;
                }
            }
        double tmpUnificationQuantity = (tmpFirstFingerprint.size()+ tmpSecondFingerprint.size()) - tmpIntersection;
        double tmpTanimotoSimilarity = tmpIntersection / tmpUnificationQuantity;
        return tmpTanimotoSimilarity;
    }
}
