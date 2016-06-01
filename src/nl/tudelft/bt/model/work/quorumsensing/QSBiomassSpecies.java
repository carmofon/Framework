package nl.tudelft.bt.model.work.quorumsensing;

import java.awt.Color;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class QSBiomassSpecies extends BiomassSpecies {
	protected SoluteSpecies _autoInducer;

	protected double _threshold;

	protected Color _colorBelowQuorum;

	protected Color _colorAboveQuorum;

	public QSBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume, SoluteSpecies autoInducer,
			double threshold, Color colorBelowQuorum, Color colorAboveQuorum)
			throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		_threshold = threshold;
		_autoInducer = autoInducer;
		_colorBelowQuorum = colorBelowQuorum;
		_colorAboveQuorum = colorAboveQuorum;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new QSBiomassParticle(this);
	}
}
