package nl.tudelft.bt.model.work.paracrineloop;

import java.awt.Color;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements the mutator particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - June 1, 2009
 */
public class InducibleBiomassParticle extends BiomassParticle {
	public InducibleBiomassParticle(InducibleBiomassSpecies s) {
		super(s);
	}

	//TODO check if these methods can be removed
	@Override
	public BiomassParticle divide() {
		InducibleBiomassParticle daughter = (InducibleBiomassParticle) super
				.divide();
		return daughter;
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	@Override
	public Color getColorCore() {
		//The problem was that setCurrentLocation was calle after getValue!
		MultigridVariable.setCurrentLocation(getCenter());
		double a = ((InducibleBiomassSpecies)_biomassSpecies)._autoInducer.getValue();
		double t = ((InducibleBiomassSpecies)_biomassSpecies)._threshold;
		Color cBelow = ((InducibleBiomassSpecies)_biomassSpecies)._colorBelowQuorum;
		Color cAbove = ((InducibleBiomassSpecies)_biomassSpecies)._colorAboveQuorum;
		
		// TODO Auto-generated method stub
		return (a > t ? cAbove : cBelow);
	}
}
