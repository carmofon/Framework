package nl.tudelft.bt.model.work.vanni;

import java.awt.Color;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.reaction.ProcessFactor;

public class BacteriocinBiomassSpecies extends BiomassSpecies {
	protected SoluteSpecies _bacteriocin;

	protected float _threshold;
	
	protected ProcessFactor _m0;

	protected Color _colorBelowQuorum;

	protected Color _colorAboveQuorum;

	public BacteriocinBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume, SoluteSpecies bacteriocin,
			float threshold, ProcessFactor m0,
			Color colorBelowQuorum, Color colorAboveQuorum)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		_threshold = threshold;
		_m0= m0;
		_bacteriocin = bacteriocin;
		_colorBelowQuorum = colorBelowQuorum;
		_colorAboveQuorum = colorAboveQuorum;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new BacteriocinBiomassParticle(this);
	}
}
