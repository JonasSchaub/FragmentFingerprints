package de.unijena.cheminf.fragment.fingerprint;

import org.openscience.cdk.fingerprint.ICountFingerprint;
import java.util.HashMap;

public class CountFingerprint implements ICountFingerprint {
    //<editor-fold desc="private  class variables" defaultstate="collapsed">
    /**
     * Is an array containing all the fragments that are to be used to generate the fingerprint.
     */
    String[] fragmentArrayOfUniqueSmiles;
    /**
     * count array
     */
    int[] countArray;
    /**
     * Map that contains all the fragments that have been removed from a molecule with their frequencies.
     */
    HashMap<String,Integer> fragmentHashMap; // TODO rename
    //</editor-fold>
    //
    /**
     * Constructor
     *
     * @param aFragments
     * @param aCountArray
     * @param aMapOfFragmentSmiles
     */
    public CountFingerprint(String[] aFragments, int[] aCountArray, HashMap<String, Integer> aMapOfFragmentSmiles) {
        this.fragmentArrayOfUniqueSmiles = aFragments;
        this.countArray = aCountArray;
        this.fragmentHashMap = aMapOfFragmentSmiles;
    }

    /**
     * Returns the size of the fingerprint, i.e., the number of hash bins.
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
    }

    /**
     * Returns the count value for the bin with the given index
     *
     * @param index the index of the bin to return the number of hits for.
     * @return int
     */
    @Override
    public int getCount(int index) {
        return  this.countArray[index];
    }

    /**
     * Returns the hash corresponding to the given index in the fingerprint
     *
     * @param index the index of the bin to return the hash for.
     * @return
     */
    @Override
    public int getHash(int index) {
        return this.fragmentHashMap.get(fragmentArrayOfUniqueSmiles[index]);
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
            return this.fragmentHashMap.containsKey(this.fragmentArrayOfUniqueSmiles[hash]);
        }
    }

    /**
     * Returns the number of times a certain hash exists in the fingerprint
     *
     * @param hash
     * @return
     */
    @Override
    public int getCountForHash(int hash) {
        return this.countArray[hash];
    }
}
