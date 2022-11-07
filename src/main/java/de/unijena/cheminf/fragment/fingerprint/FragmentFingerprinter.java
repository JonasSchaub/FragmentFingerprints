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

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Class to generate a bit fragment fingerprint
 *
 * @author Betuel Sevindik
 */
public class FragmentFingerprinter implements IFragmentFingerprinter {

    BitSet bitSet;
    int[] bitVector;

    int[] countVector;
    ArrayList<String> fragmentSetList;
    //</editor-fold>
    //
    /**
     * Constructor
     */
    public FragmentFingerprinter(ArrayList<String> aFragmentList) {
        this.fragmentSetList = aFragmentList;
    }

    /**
     * Empty Constructor
     */
    public FragmentFingerprinter() {

    }
    //
    //<editor-fold desc="public methods" defaultstate="collapsed">
    /**
     *
     *
     * @param aListWithUniqueSmiles
     * @return
     */
    @Override
    public IBitFingerprint getBitFingerprint(ArrayList<String> aListWithUniqueSmiles) {
        int tmpVectorSize = this.fragmentSetList.size();
        this.bitVector = new int[tmpVectorSize];
        this.bitSet = new BitSet(this.fragmentSetList.size()); // TODO
        for (int tmpDefaultBit = 0; tmpDefaultBit < tmpVectorSize; tmpDefaultBit++) {
            this.bitVector[tmpDefaultBit] = 0;
        }
        /**
        String[] keys = new String[aMoleculeFragmentsMap.size()];
        int iter = 0;
        for(String tmpKey : aMoleculeFragmentsMap.keySet()) {
            keys[iter]  = tmpKey;
            iter++;
        }
         */
        HashMap<String, Integer> tmpFragmentHashMap = new HashMap<>();
        // sort aFragmentList alphabetically
        Collections.sort(this.fragmentSetList, String.CASE_INSENSITIVE_ORDER);
        int tmpValuePosition = 0;
        for (String tmpKey : this.fragmentSetList) {
            tmpFragmentHashMap.put(tmpKey, tmpValuePosition);
            tmpValuePosition++;
        }
        for (String tmpUniqueSmiles : aListWithUniqueSmiles) {
            if (tmpFragmentHashMap.containsKey(tmpUniqueSmiles)) {
                int tmpPosition = tmpFragmentHashMap.get(tmpUniqueSmiles);
                this.bitVector[tmpPosition] = 1;
                this.bitSet.set(tmpPosition,true);
            }
        }
        return new IBitFingerprint() {
            @Override
            public int cardinality() {
                return bitSet.cardinality();
            }
            @Override
            public long size() {
                return bitSet.size();
            }
            @Override
            public void and(IBitFingerprint fingerprint) {
                throw new UnsupportedOperationException();
            }
            @Override
            public void or(IBitFingerprint fingerprint) {
                throw new UnsupportedOperationException();
            }
            @Override
            public boolean get(int index) {
                return bitSet.get(index);
            }
            @Override
            public void set(int index, boolean value) {
                throw new UnsupportedOperationException();
            }
            @Override
            public BitSet asBitSet() {
                return (BitSet) bitSet.clone();
            }
            @Override
            public void set(int i) {
                throw new UnsupportedOperationException();
            }
            @Override
            public int[] getSetbits() {
                int[] tmpPositiveBits = new int[bitSet.cardinality()];
                int tmpIndex = 0;
                for (int i = 0; i < bitSet.length(); i++) {
                    if (bitSet.get(i)) {
                        tmpPositiveBits[tmpIndex++] = i;
                    }
                }
                return tmpPositiveBits;
            }
        };
    }
    @Override
    public ICountFingerprint getCountFingerprint(HashMap<String, Integer> aMoleculeFragmentsMap) {
        int tmpVectorSize = this.fragmentSetList.size();
        this.countVector = new int[tmpVectorSize];
        // Create null vector
        for (int tmpDefaultBit = 0; tmpDefaultBit < tmpVectorSize; tmpDefaultBit++) {
            this.countVector[tmpDefaultBit] = 0;
        }
        HashMap<String, Integer> tmpFragmentHashMap = new HashMap<>();
        // Sort aFragmentList alphabetically
        Collections.sort(this.fragmentSetList, String.CASE_INSENSITIVE_ORDER);
        int tmpValuePosition = 0;
        // Convert the fragment list in a HashMap
        for (String tmpKey : this.fragmentSetList) {
            tmpFragmentHashMap.put(tmpKey, tmpValuePosition);
            tmpValuePosition++;
        }
        /*
         *  Comparison of two hashmaps on identical keys
         *  TODO: Possibly make it more efficient
         */
        Map<Integer, Integer> tmpCountMap = new HashMap<>(this.fragmentSetList.size());
        for (String tmpUniqueSmiles : tmpFragmentHashMap.keySet()) {
            if (aMoleculeFragmentsMap.containsKey(tmpUniqueSmiles)) {
                int tmpPosition = tmpFragmentHashMap.get(tmpUniqueSmiles);
                this.countVector[tmpPosition] = aMoleculeFragmentsMap.get(tmpUniqueSmiles);
                //set.set(Integer.parseInt(aMoleculeFragments.get(tmpUniqueSmiles)));
                tmpCountMap.put(tmpFragmentHashMap.get(tmpUniqueSmiles),aMoleculeFragmentsMap.get(tmpUniqueSmiles));
            }
            else{
                tmpCountMap.put(tmpFragmentHashMap.get(tmpUniqueSmiles), 0);
            }
        }
        int[] tmpHash = new int[tmpCountMap.size()];
        int[] tmpCount = new int[tmpCountMap.size()];
        int i = 0;
        for(int h: tmpCountMap.keySet()) {
            tmpHash[i] = h;
            tmpCount[i++] = tmpCountMap.get(h);
        }
        return new ICountFingerprint() {
            @Override
            public long size() {
                return fragmentSetList.size();
            }
            @Override
            public int numOfPopulatedbins() {
                return tmpCountMap.size();
            }
            @Override
            public int getCount(int index) {
                return tmpCount[index];
            }
            @Override
            public int getHash(int index) {
                return tmpHash[index];
            }
            @Override
            public void merge(ICountFingerprint fp) {
                throw new UnsupportedOperationException();
            }
            @Override
            public void setBehaveAsBitFingerprint(boolean behaveAsBitFingerprint) {
                throw new UnsupportedOperationException();
            }
            @Override
            public boolean hasHash(int hash) {
                return tmpCountMap.containsKey(hash);
            }
            @Override
            public int getCountForHash(int hash) {
                return tmpCountMap.get(hash);
            }
        };
    }

    @Override
    public String getVersionDescription() {
       throw new UnsupportedOperationException();
    }

    @Override
    public BitSet getFingerprint(IAtomContainer mol) throws CDKException {
       throw new UnsupportedOperationException();
    }

    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
        throw new UnsupportedOperationException();
    }

    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer container) throws CDKException {
        throw new UnsupportedOperationException();
    }

    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer container) throws CDKException {
       throw new UnsupportedOperationException();
    }

    @Override
    public int getSize() {
     return this.fragmentSetList.size();
    }
    public int[] getCountVector() {
        return this.countVector;
    }
    public int[] getBitVector() {
        return this.bitVector;
    }
}
