package nl.tudelft.bt.model.multigrid.boundary_layers;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.InvalidValueException;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;
/**
 * Implements moving plannar boundary layer (set at a guiven height from the
 * heighest biofilm feature)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class FrontBoundaryLayer extends BoundaryLayer {
	// height above the biofilm where boundary layer is placed
	private float _maximumH;
	private float _h;
	/**
	 * Constructs a new boundary layer object with planar boundary layer placed
	 * h um above the biofilm surface
	 * 
	 * @param h
	 */
	public FrontBoundaryLayer(float h) throws ModelException {
		super();
		// check validity of h parameter, i.e. maximum biofilm height
		// plus h must be less than the system height.
		if ((Model.model().getMaximumBiofilmHeight() + h) >= Model.model().systemSize.x)
			throw new InvalidValueException(
					"H value for front boundary layer (" + h + ")is to high");
		_h = h;
		_maximumH = h;
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.multigrid.BoundaryLayer#setBoundaryLayer(org.photobiofilms.phlip.multigrid.BacteriaSpecies[])
	 */
	public void setBoundaryLayer(ParticulateSpecies[] b,
			BoundaryConditions bc) {
		float[][][] bl = _mg[_order - 1];
		int n = bl.length;
		int m = bl[0].length;
		int l = bl[0][0].length;
		// find the boundary layer height
		int h = (int) ((Model.model().getCurrentBiofilmHeight() + _h) / _voxelSide) + 1;
		for (int i = 0; i < n; i++) {
			float v = (i > h ? 1.0f : 0.0f);
			for (int j = 0; j < m; j++)
				for (int k = 0; k < l; k++)
					bl[i][j][k] = v;
		}
	}
	/**
	 * @param h
	 *            The _h to set.
	 */
	public void setThickness(float h) {
		//		if (h > _maximumH) {
		//			throw new ModelRuntimeException("Trying to set heigh of a"
		//					+ " FrontBoundaryLayer to an illegal value");
		//		} else {
		this._h = h;
		//		}
	}
}