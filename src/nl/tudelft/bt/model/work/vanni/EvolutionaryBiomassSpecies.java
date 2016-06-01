package nl.tudelft.bt.model.work.vanni;

import java.awt.Color;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class EvolutionaryBiomassSpecies extends BiomassSpecies {
	private ParticulateSpecies _parameter;
	protected float _mutationRate;
	protected float _growthfactor;
	private Color _color;

	public EvolutionaryBiomassSpecies(String name, ParticulateSpecies[] species,
			float[] fractionalCompositionInVolume, float mutationRate)
	throws NonMatchingNumberException {
		super(name, species, fractionalCompositionInVolume);
		_mutationRate = mutationRate;
	}
	
	public void setParameterValue(ParticulateSpecies s)
			throws ModelRuntimeException {
		for (int i = 0; i < _particulates.length; i++) {
			if (s == _particulates[i]) {
				_parameter = s;
				return;
			}
		}
		throw new ModelRuntimeException("Attempt to set non-existing species "
				+ s.getName() + " as active mass in " + getName());
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassSpecies#createBiomassParticle()
	 */
	public BiomassParticle createBiomassParticle() {
		return new EvolutionaryBiomassParticle(this);
	}

	public float getGrowthFactor() {
		// TODO Auto-generated method stub
		_growthfactor = 1;
		return _growthfactor;
	}

}
