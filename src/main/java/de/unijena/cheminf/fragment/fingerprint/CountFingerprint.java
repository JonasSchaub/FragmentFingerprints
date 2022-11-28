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

import org.openscience.cdk.fingerprint.ICountFingerprint;

import java.util.HashMap;

/**
 * The CountFingerprint class implements the CDK interface ICountFingerprint.
 * ICountFingerprint provides useful methods to obtain information about the calculated count fingerprint.
 */
public class CountFingerprint implements ICountFingerprint {
    //<editor-fold desc="private  class variables" defaultstate="collapsed">
    /**
     * Is an array containing all predefined unique SMILES.
     */
    String[] fragmentArrayOfUniqueSmiles;
    /**
     * The HashMap maps the unique SMILES to the position they have in the array.
     */
    HashMap<String,Integer> uniqueSmilesToPositionMap;
    /**
     * Map
     */
    HashMap<Integer, Integer> rawMap;
    //</editor-fold>
    //
    /**
     * Constructor.
     *
     *
     * @param aFragments
     * @param aMapOfFragmentSmiles
     */
    public CountFingerprint(String[] aFragments, HashMap<String, Integer> aMapOfFragmentSmiles, HashMap<Integer, Integer> aRawMap) {
        this.fragmentArrayOfUniqueSmiles = aFragments;
        this.uniqueSmilesToPositionMap = aMapOfFragmentSmiles;
        this.rawMap = aRawMap;
    }

    /**
     * Returns the number of bits in this fingerprint.
     * Since fragment fingerprints are key based, the number of bits in the fingerprint
     * is equal to the number of predefined fragments (unique SMILES).
     *
     * @return the size of the fingerprint.
     */
    @Override
    public long size() {
        return this.fragmentArrayOfUniqueSmiles.length;
    }

    /**
     * Returns the number of bins that are populated.
     * The number is equal to the defined fragments.
     *
     * @return the number of populated bins
     * @see    #size()
     */
    @Override
    public int numOfPopulatedbins() {
        return this.fragmentArrayOfUniqueSmiles.length;
    }

    /**
     * For each defined structure in the fingerprint, returns the frequency.
     * For an index >= 0 and < fragmentArrayOfUniqueSmiles, the frequency of this unique SMILES is output.
     * If the given index is > than the size of the fingerprint, a NullPinterException is thrown.
     *
     * @param index the index of the bin to return the number of hits for.
     * @return the count for the bin with given index.
     */
    @Override
    public int getCount(int index) {
        if (index < this.fragmentArrayOfUniqueSmiles.length && index >= 0 && this.rawMap.containsKey(index)) {
            return this.rawMap.get(index);
        } else if (index >= this.fragmentArrayOfUniqueSmiles.length) {
            return this.rawMap.get(index);
        } else {
            return 0;
        }
    }

    /**
     * Returns the hash corresponding to the given index in the fingerprint.
     * Since this is a key-based fingerprint, the hash value is nothing
     * more than the position of the bin in the fingerprint (no hash value is calculated).
     *
     * @param index the index of the bin to return the hash for.
     * @return the hash for the bin with the given index.
     */
    @Override
    public int getHash(int index) {
        return this.uniqueSmilesToPositionMap.get(this.fragmentArrayOfUniqueSmiles[index]);
    }

    /**
     *
     * @param fp to be merged
     */
    @Override
    public void merge(ICountFingerprint fp) {
        throw new UnsupportedOperationException();
    }

    /**
     *
     * @param behaveAsBitFingerprint
     */
    @Override
    public void setBehaveAsBitFingerprint(boolean behaveAsBitFingerprint) {
        throw new UnsupportedOperationException();
    }

    /**
     * Whether the fingerprint contains the given hash.
     *The parameter hash is not a calculated hash value, but also corresponds to the position of the bin in the fingerprint.
     * @see    #getHash(int)
     * For hash >= fingerprint size the method returns false and for all hash < fingerprint size a true is returned.
     *
     * @param hash
     * @return true if the fingerprint contains the given hash, otherwise false.
     */
    @Override
    public boolean hasHash(int hash) {
        if(hash >= this.fragmentArrayOfUniqueSmiles.length){
            return false;
        } else {
            return this.uniqueSmilesToPositionMap.containsKey(this.fragmentArrayOfUniqueSmiles[hash]);
        }
    }

    /**
     * Get the number of times a certain hash exists in the fingerprint.
     * Since the fragment fingerprint is a key-based fingerprint and the hash value therefore indicates the
     * position of the respective bin in the fingerprint, each hash value occurs only once in the fingerprint.
     * Therefore, the method returns the frequency of the fragment in the specified bin.
     *
     * @param hash the index of the bin to return the number of hits for.
     * @return the number associated with the given hash
     */
    @Override
    public int getCountForHash(int hash) {
        if (hash < this.fragmentArrayOfUniqueSmiles.length && hash >= 0 && this.rawMap.containsKey(hash)) {
            return this.rawMap.get(hash);
        } else if (hash >= this.fragmentArrayOfUniqueSmiles.length) {
            return this.rawMap.get(hash);
        } else {
            return 0;
        }
    }
}
