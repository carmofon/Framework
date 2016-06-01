package nl.tudelft.bt.model.particlebased.granule.phbgranule;

import java.io.Serializable;

import org.jfree.chart.ui.ColorPalette;

/**
 * A color palette with gray tones.
 *
 * @author Joao Xavier
 */
public class BonePalette extends ColorPalette implements Serializable {

    /**
     * Creates a new palette.
     */
    public BonePalette() {
        super();
        initialize();
    }

    /**
     * Intializes the palette's indices.
     */
    public void initialize() {

        setPaletteName("Bone");

        r = new int[256];
        g = new int[256];
        b = new int[256];


        for (int i = 0; i < 256; i++) {
        		int v = (int)(((float)i)/255f * 200) + 20;
            r[i] = v;
            g[i] = v;
            b[i] = v;
        }

    }

}
