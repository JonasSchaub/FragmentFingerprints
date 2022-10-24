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
import java.util.HashMap;

/**
 * Interface for implementing fragment fingerprints
 *
 * @author Betuel Sevindik
 */
public interface IFingerprint {

    // <editor-fold defaultstate="collapsed" desc="Methods">
    /**
     * Returns a List with the fingerprint and the positive indices in the fingerprint array
     *
     * @param aList  a fragment list
     * @param aMap  a molecule HashMap
     * @return list
     */
    ArrayList generateFragmentFingerprint(ArrayList<String> aList, HashMap<String,String> aMap);
    //
    /**
     * Returns the fingerprint size
     *
     * @return int number
     */
    int getFragmentFingerprintSize();
    //
    /**
     * Returns the value of the similarity between two fingerprints
     *
     * @return
     */
    double calculateTanimotoSimilarity();
    //
    /**
     * Returns the number of positive bits in the fingerprint
     *
     * @return
     */
    int getNumberPositiveBits();
    //
    /**
     * Returns the positive indices
     *
     * @param aList
     * @return
     */
    Object getIndicesPositiveBits(ArrayList<Object> aList);
    // </editor-fold>
    //

    //TODO sparse vector

















}
