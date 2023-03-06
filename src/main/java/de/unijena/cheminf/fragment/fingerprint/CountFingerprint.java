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
 */

package de.unijena.cheminf.fragment.fingerprint;

import org.openscience.cdk.fingerprint.ICountFingerprint;
import java.util.HashMap;

/**
 * The CountFingerprint class implements the CDK interface ICountFingerprint.
 * ICountFingerprint provides useful methods to obtain information about the calculated count fingerprint.
 *
 * @author Betuel Sevindik
 */
public class CountFingerprint implements ICountFingerprint {
    //<editor-fold desc="private final class variables" defaultstate="collapsed">
    /**
     * An array containing all predefined unique SMILES.
     */
    private final String[] fragmentArrayOfUniqueSmiles;
    /**
     * The HashMap maps the unique SMILES to the position they have in the array.
     */
    private final HashMap<String,Integer> uniqueSmilesToPositionMap;
    /**
     * The HashMap maps the position of unique SMILES to the frequency
     */
    private final HashMap<Integer, Integer> uniqueSmilesPositionToFrequencyCountRawMap;
    //</editor-fold>
    //
    //<editor-fold desc="private class variables" defaultstate="collapsed">
    private boolean behaveAsBitFingerprint;
    //</editor-fold>
    //
    //<editor-fold desc="Constructor" defaultstate="collapsed">
    /**
     * Constructor.
     *
     * @param aFragments is a string array that stores all fragments that are in the form of unique SMILES.
     * @param aSmilesToFragmentMap is a map that assigns the unique SMILES to the positions they occupy in aFragments.
     * @param aPositionToFrequencyMap is the raw map. It maps the positions of the unique SMILES to the
     *                                frequency of these SMILES.
     */
    public CountFingerprint(String[] aFragments, HashMap<String, Integer> aSmilesToFragmentMap, HashMap<Integer, Integer> aPositionToFrequencyMap) {
        if(aFragments == null || aSmilesToFragmentMap == null || aPositionToFrequencyMap == null) {
            throw new NullPointerException("At least one of the arguments is null.");
        }
        this.fragmentArrayOfUniqueSmiles = aFragments;
        this.uniqueSmilesToPositionMap = aSmilesToFragmentMap;
        this.uniqueSmilesPositionToFrequencyCountRawMap = aPositionToFrequencyMap;
        this.behaveAsBitFingerprint = false;
    }
    //</editor-fold>
    //
    //<editor-fold desc="Overridden public methods " defaultstate="collapsed">
    /**
     * {@inheritDoc}
     *
     * Since fragment fingerprints are key-based, the number of bits in the fingerprint
     * is equal to the number of predefined fragments (unique SMILES).
     *
     */
    @Override
    public long size() {
        return this.fragmentArrayOfUniqueSmiles.length;
    }
    //
    /**
     * {@inheritDoc}
     *
     * Fragment fingerprints are key-based fingerprints,
     * therefore the number of populated bins corresponds to the number of predefined fragments (unique SMILES).
     *
     */
    @Override
    public int numOfPopulatedbins() {
        return this.fragmentArrayOfUniqueSmiles.length;
    }
    //
    /**
     * {@inheritDoc}
     *
     * For an index specified within the fingerprint size, the corresponding frequency is output.
     * If the given index is greater than the size of the fingerprint or a negative value,
     * an IllegalArgumentException is thrown.
     *
     * @throws IllegalArgumentException is thrown if the given index do not exist in the fingerprint.
     */
    @Override
    public int getCount(int index) throws IllegalArgumentException {
        if(this.behaveAsBitFingerprint == false) {
            if (index >= 0 && this.uniqueSmilesPositionToFrequencyCountRawMap.containsKey(index)) {
                return this.uniqueSmilesPositionToFrequencyCountRawMap.get(index);
            } else if (index >= this.fragmentArrayOfUniqueSmiles.length || index < 0) {
                throw new IllegalArgumentException("This position does not exist in the fingerprint (undefined state).");
            } else {
                return 0;
            }
        } else {
            if (index >= 0 && this.uniqueSmilesPositionToFrequencyCountRawMap.containsKey(index)) {
                return 1;
            } else if (index >= this.fragmentArrayOfUniqueSmiles.length || index < 0) {
                throw new IllegalArgumentException("This position does not exist in the fingerprint ( undefined state).");
            } else {
                return 0;
            }
        }
    }
    //
    /**
     * {@inheritDoc}
     *
     * Since this is a key-based fingerprint, the hash value is  the position of the bin in the
     * fingerprint (no hash value is calculated).
     *
     */
    @Override
    public int getHash(int index) {
        return this.uniqueSmilesToPositionMap.get(this.fragmentArrayOfUniqueSmiles[index]);
    }
    //
    /**
     * {@inheritDoc}
     */
    @Override
    public void merge(ICountFingerprint fp) {
        throw new UnsupportedOperationException();
    }
    //
    /**
     * {@inheritDoc}
     */
    @Override
    public void setBehaveAsBitFingerprint(boolean behaveAsBitFingerprint) {
        this.behaveAsBitFingerprint = behaveAsBitFingerprint;
    }
    //
    /**
     * {@inheritDoc}
     *
     * The parameter hash is not a calculated hash value, but also corresponds to the
     * position of the bin in the fingerprint.
     *
     * @throws IllegalArgumentException is thrown if the given hash value is negative.
     */
    @Override
    public boolean hasHash(int hash) throws IllegalArgumentException {
        if(hash< this.fragmentArrayOfUniqueSmiles.length && hash>=0) {
            return true;
        } else if (hash < 0) {
            throw new IllegalArgumentException("Negative values are not allowed.");
        } else {
            return false;
        }
    }
    //
    /**
     * {@inheritDoc}
     *
     * Since the fragment fingerprint is a key-based fingerprint and the hash value therefore indicates the
     * position of the respective bin in the fingerprint, each hash value occurs only once in the fingerprint.
     * Therefore, the method returns the frequency of the fragment in the specified bin.
     *
     * @see #getCount(int)
     *
     * @throws IllegalArgumentException is thrown if the given hash value do not exist in the fingerprint.
     *
     */
    @Override
    public int getCountForHash(int hash) throws IllegalArgumentException {
        if(this.behaveAsBitFingerprint == false) {
            if (hash >= 0 && this.uniqueSmilesPositionToFrequencyCountRawMap.containsKey(hash)) {
                return this.uniqueSmilesPositionToFrequencyCountRawMap.get(hash);
            } else if (hash >= this.fragmentArrayOfUniqueSmiles.length || hash < 0) {
                throw new IllegalArgumentException("This position does not exist in the fingerprint (undefined state).");
            } else {
                return 0;
            }
        } else {
            if (hash >= 0 && this.uniqueSmilesPositionToFrequencyCountRawMap.containsKey(hash)) {
                return 1;
            } else if (hash >= this.fragmentArrayOfUniqueSmiles.length || hash < 0) {
                throw new IllegalArgumentException("This position does not exist in the fingerprint (undefined state).");
            } else {
                return 0;
            }
        }
    }
    //</editor-fold>
    //
}
