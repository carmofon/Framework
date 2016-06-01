package nl.tudelft.bt.model.work.vanni;

import java.awt.Color;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

public class EvolutionaryBiomassParticle extends BiomassParticle {
	protected float _growthfactor;
	private EvolutionaryBiomassSpecies _species;
    public Color _color;

	public EvolutionaryBiomassParticle(EvolutionaryBiomassSpecies s) {
		super(s);
		// TODO pseudocode
		_growthfactor = s.getGrowthFactor();
		_species=s;
	}

	public EvolutionaryBiomassParticle(EvolutionaryBiomassSpecies s, float P) {
		super(s);
		_growthfactor = P * Model.model().getRandomFromNormalDistribution();
	}

	@Override
	public EvolutionaryBiomassParticle divide() {
		EvolutionaryBiomassParticle daughter = (EvolutionaryBiomassParticle) super.divide();
		if (_species._mutationRate > Model.model().getRandom()) {
			//assign a new growth factor randomly drawn
			float newGrowthFactor = _growthfactor * (Model.model().getRandomFromNormalDistribution() + 1);
			daughter.setGrowthFactor(newGrowthFactor);
			daughter.getColorCore();
		}else{
			float newGrowthFactor = _growthfactor * 1;	
			daughter.setGrowthFactor(newGrowthFactor);
			daughter.getColorCore();
		}

		return daughter;
	}

	protected void setGrowthFactor(float P) {
		_growthfactor = P;
	}

	@Override
	public float grow(float t){ 
		float m = super._composition.grow(t, getCenter());
		m = m * this._growthfactor;
		return m;
	}
	
	public Color getColorCore() {
		
		//float iRed = ((float) (_inducedColor.getRed())) / 255f;
		//float iGreen = ((float) (_inducedColor.getGreen())) / 255f;
		//float iBlue = ((float) (_inducedColor.getBlue())) / 255f;
		//float ColorArray[];
		float iRed=(1.0f/(1.0f+this._growthfactor));
		float iGreen=(1.0f/(1.0f+this._growthfactor));
		float iBlue=(1.0f/(1.0f+this._growthfactor));
		return new Color(iRed , iGreen, iBlue);
	}
}
