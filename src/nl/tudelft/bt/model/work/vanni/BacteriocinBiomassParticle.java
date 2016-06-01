package nl.tudelft.bt.model.work.vanni;

import java.awt.Color;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.work.quorumsensing.QSBiomassSpecies;

/**
 * Implements the mutator particle
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 18, 2006
 */
public class BacteriocinBiomassParticle extends BiomassParticle {
	public BacteriocinBiomassParticle(BacteriocinBiomassSpecies s) {
		super(s);
	}

	//TODO check if these methods can be removed
	@Override
	public BiomassParticle divide() {
		BacteriocinBiomassParticle daughter = (BacteriocinBiomassParticle) super
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
		double a = ((BacteriocinBiomassSpecies)_biomassSpecies)._bacteriocin.getValue();
		double t = ((BacteriocinBiomassSpecies)_biomassSpecies)._threshold;
		t = t * ((BacteriocinBiomassSpecies)_biomassSpecies)._m0.getValue();
		Color cBelow = ((BacteriocinBiomassSpecies)_biomassSpecies)._colorBelowQuorum;
		Color cAbove = ((BacteriocinBiomassSpecies)_biomassSpecies)._colorAboveQuorum;
		
		// TODO Auto-generated method stub
		return (a > t ? cAbove : cBelow);
	}
}
