package de.unijena.cheminf.fragment.fingerprint;

import org.openscience.cdk.fingerprint.ICountFingerprint;

import java.util.HashMap;
import java.util.Map;

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
     * Returns the number of bits of this fingerprint.
     *
     * @return
     */
    @Override
    public long size() {
        return this.fragmentArrayOfUniqueSmiles.length;
    }

    /**
     * Returns the number of bins that are populated
     *
     * @return int
     */
    @Override
    public int numOfPopulatedbins() {
        return this.fragmentArrayOfUniqueSmiles.length;
    } // TODO

    /**
     * Returns the count value for the bin with the given index
     *
     * @param index the index of the bin to return the number of hits for.
     * @return int
     */
    @Override
    public int getCount(int index)  {
            if (index < this.fragmentArrayOfUniqueSmiles.length && this.rawMap.containsKey(index)) { // TODO größer null
                return this.rawMap.get(index);
            } else if (index >= this.fragmentArrayOfUniqueSmiles.length) {
                return this.rawMap.get(index);
            } else {
                return 0;
            }
    } //this.countArray[index]

    /**
     * Returns the hash corresponding to the given index in the fingerprint
     *
     * @param index the index of the bin to return the hash for.
     * @return
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
     * Whether the fingerprint contains the given hash
     *
     * @param hash
     * @return  boolean (true or false)
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
     * Returns the number of times a certain hash exists in the fingerprint
     *
     * @param hash
     * @return
     */
    @Override
    public int getCountForHash(int hash) {  //return this.countArray[hash];
      throw new UnsupportedOperationException();
    }
}
