/*
 * @(#)AppleRLEEncoder.java  1.0  2011-01-05
 * 
 * Copyright © 2011 Werner Randelshofer, Immensee, Switzerland.
 * All rights reserved.
 * 
 * You may not use, copy or modify this file, except in compliance with the
 * license agreement you entered into with Werner Randelshofer.
 * For details see accompanying license terms.
 */
package ch.randelshofer.media.quicktime;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.imageio.stream.ImageOutputStream;
import javax.imageio.stream.ImageOutputStreamImpl;
import javax.imageio.stream.MemoryCacheImageOutputStream;

/**
 * Implements the run length encoding of the Apple QuickTime Animation (RLE)
 * format.
 * <p>
 * An RLE-encoded frame has the following format:
 * <p>
 * <pre>
 * Header:
 * uint32 chunkSize
 *
 * uint16 header 0x0000 => decode entire image
 *               0x0008 => starting line and number of lines follows
 * if header==0x0008 {
 *   uint16 startingLine at which to begin updating frame
 *   uint16 reserved 0x0000
 *   uint16 numberOfLines to update
 *   uint16 reserved 0x0000
 * }
 * n-bytes compressed lines
 * </pre>
 *
 * The first 4 bytes defines the chunk length. This field also carries some
 * other unknown flags, since at least one of the high bits is sometimes set.<br>
 *
 * If the overall length of the chunk is less than 8, treat the frame as a
 * NOP, which means that the frame is the same as the one before it.<br>
 *
 * Next, there is a header of either 0x0000 or 0x0008. A header value with
 * bit 3 set (header &amp; 0x0008) indicates that information follows revealing
 * at which line the decode process is to begin:<br>
 *
 * <pre>
 * 2 bytes    starting line at which to begin updating frame
 * 2 bytes    unknown
 * 2 bytes    the number of lines to update
 * 2 bytes    unknown
 * </pre>
 *
 * If the header is 0x0000, then the decode begins from the first line and
 * continues through the entire height of the image.<br>
 *
 * After the header comes the individual RLE-compressed lines. An individual
 * compressed line is comprised of a skip code, followed by a series of RLE
 * codes and pixel data:<br>
 * <pre>
 *  1 byte     skip code
 *  1 byte     RLE code
 *  n bytes    pixel data
 *  1 byte     RLE code
 *  n bytes    pixel data
 * </pre>
 * After the skip byte is the first RLE code, which is a single signed
 * byte. The RLE code can have the following meanings:<br>
 * <ul>
 * <li>equal to 0: There is another single-byte skip code in the stream.
 *              Again, the actual number of pixels to skip is 1 less
 *              than the skip code. Therefore, the maximum skip byte
 *              value of 255 allows for a maximum of 254 pixels to be
 *              skipped.</li>
 *
 * <li>equal to -1: End of the RLE-compressed line</li>
 *
 * <li>greater than 0: Run of pixel data is copied directly from the
 *              encoded stream to the output frame.</li>
 *
 * <li>less than -1: Repeat pixel data -(RLE code) times.</li>
 * </ul>
 * <p>
 * The pixel data has the following format:
 * <ul>
 * <li>8-bit data: Pixels are handled in groups of four. Each pixel is a palette
 * index (the palette is determined by the Quicktime file transporting the
 * data).<br>
 * If (code &gt; 0), copy (4 * code) pixels from the encoded stream to the
 * output.<br>
 * If (code &lt; -1), extract the next 4 pixels from the encoded stream
 * and render the entire group -(code) times to the output frame. </li>
 *
 * <li>16-bit data: Each pixel is represented by a 16-bit RGB value with 5 bits
 * used for each of the red, green, and blue color components and 1 unused bit
 * to round the value out to 16 bits: {@code xrrrrrgg gggbbbbb}. Pixel data is
 * rendered to the output frame one pixel at a time.<br>
 * If (code &gt; 0), copy the run of (code) pixels from the encoded stream to
 * the output.<br>
 * If (code &lt; -1), unpack the next 16-bit RGB value from the encoded stream
 * and render it to the output frame -(code) times.</li>
 *
 * <li>24-bit data: Each pixel is represented by a 24-bit RGB value with 8 bits
 * (1 byte) used for each of the red, green, and blue color components:
 * {@code rrrrrrrr gggggggg bbbbbbbb}. Pixel data is rendered to the output
 * frame one pixel at a time.<br>
 * If (code &gt; 0), copy the run of (code) pixels from the encoded stream to
 * the output.<br>
 * If (code &lt; -1), unpack the next 24-bit RGB value from the encoded stream
 * and render it to the output frame -(code) times.</li>
 *
 * <li>32-bit data: Each pixel is represented by a 32-bit ARGB value with 8 bits
 * (1 byte) used for each of the alpha, red, green, and blue color components:
 * {@code aaaaaaaa rrrrrrrr gggggggg bbbbbbbb}. Pixel data is rendered to the
 * output frame one pixel at a time.<br>
 * If (code &gt; 0), copy the run of (code) pixels from the encoded stream to
 * the output.<br>
 * If (code &lt; -1), unpack the next 32-bit ARGB value from the encoded stream
 * and render it to the output frame -(code) times.</li>
 * </ul>
 *
 * References:<br/>
 * <a href="http://multimedia.cx/qtrle.txt">http://multimedia.cx/qtrle.txt</a><br>
 *
 * @author Werner Randelshofer
 * @version 1.0 2011-01-05 Created.
 */
public class AppleRLEEncoder {

    /** Encodes a 16-bit key frame.
     * 
     * @param out The output stream. Must be set to Big-Endian.
     * @param data The image data.
     * @param offset The offset to the first pixel in the data array.
     * @param length The width of the image in data elements.
     * @param step The number to add to offset to get to the next scanline.
     */
    public void writeKey16(ImageOutputStream out, short[] data, int offset, int length, int step)
            throws IOException {
        long headerPos = out.getStreamPosition();

        // Reserve space for the header:
        out.writeInt(0);
        out.writeShort(0x0000);

        // Encode each scanline
        int ymax = data.length;
        for (int y = offset; y < ymax; y += step) {
            int xy = y;
            int xymax = y + length;

            out.write(1); // this is a key-frame, there is nothing to skip at the start of line

            int literalCount = 0;
            int repeatCount = 0;
            for (; xy < xymax; ++xy) {
                // determine repeat count 
                short v = data[xy];
                for (repeatCount = 0; xy < xymax && repeatCount < 127; ++xy, ++repeatCount) {
                    if (data[xy] != v) {
                        break;
                    }
                }
                xy -= repeatCount;

                if (repeatCount < 3) {
                    literalCount++;
                    if (literalCount > 126) {
                        out.write(literalCount); // Literal OP-code
                        out.writeShorts(data, xy - literalCount + 1, literalCount);
                        literalCount = 0;
                    }
                } else {
                    if (literalCount > 0) {
                        out.write(literalCount); // Literal OP-code
                        out.writeShorts(data, xy - literalCount, literalCount);
                        literalCount = 0;
                    }
                    out.write(-repeatCount); // Repeat OP-code
                    out.writeShort(v);
                    xy += repeatCount - 1;
                }
            }

            // flush literal run
            if (literalCount > 0) {
                out.write(literalCount);
                out.writeShorts(data, xy - literalCount, literalCount);
                literalCount = 0;
            }

            out.write(-1);// End of line OP-code
        }


        // Complete the header
        long pos = out.getStreamPosition();
        out.seek(headerPos);
        out.writeInt((int) (pos - headerPos));
        out.seek(pos);
    }

    /** Encodes a 16-bit delta frame.
     *
     * @param out The output stream. Must be set to Big-Endian.
     * @param data The image data.
     * @param prev The image data of the previous frame.
     * @param offset The offset to the first pixel in the data array.
     * @param length The width of the image in data elements.
     * @param step The number to add to offset to get to the next scanline.
     */
    public void writeDelta16(ImageOutputStream out, short[] data, short[] prev, int offset, int length, int step)
            throws IOException {

        int height = (data.length - offset) / step;

        // Determine whether we can skip lines at the beginning
        int ymin;
        int ymax = offset + height * step;
        scanline:
        for (ymin = offset; ymin < ymax; ymin += step) {
            int xy = ymin;
            int xymax = ymin + length;
            for (; xy < xymax; ++xy) {
                if (data[xy] != prev[xy]) {
                    break scanline;
                }
            }
        }


        if (ymin == ymax) {
            // => Frame is identical to previous one
            out.writeInt(4);
            return;
        }

        // Determine whether we can skip lines at the end
        scanline:
        for (; ymax > ymin; ymax -= step) {
            int xy = ymax - step;
            int xymax = ymax - step + length;
            for (; xy < xymax; ++xy) {
                if (data[xy] != prev[xy]) {
                    break scanline;
                }
            }
        }
        //System.out.println("AppleRLEEncoder ymin:" + ymin / step + " ymax" + ymax / step);

        // Reserve space for the header
        long headerPos = out.getStreamPosition();
        out.writeInt(0);

        if (ymin == offset && ymax == offset + height * step) {
            // => we can't skip any lines:
            out.writeShort(0x0000);
        } else {
            // => we can skip lines:
            out.writeShort(0x0008);
            out.writeShort(ymin / step);
            out.writeShort(0);
            out.writeShort((ymax - ymin + 1) / step);
            out.writeShort(0);
        }


        // Encode each scanline
        for (int y = ymin; y < ymax; y += step) {
            int xy = y;
            int xymax = y + length;

            // determine skip count
            int skipCount = 0;
            for (; xy < xymax; ++xy, ++skipCount) {
                if (data[xy] != prev[xy]) {
                    break;
                }
            }
            if (skipCount == length) {
                // => the entire line can be skipped
                out.write(1); // don't skip any pixels
                out.write(-1); // end of line
                continue;
            }
            out.write(Math.min(255, skipCount + 1));
            if (skipCount > 254) {
                skipCount -= 254;
                while (skipCount > 254) {
                    out.write(0); // Skip OP-code
                    out.write(255);
                    skipCount -= 254;
                }
                out.write(0); // Skip OP-code
                out.write(skipCount + 1);
            }

            int literalCount = 0;
            int repeatCount = 0;
            for (; xy < xymax; ++xy) {
                // determine skip count
                for (skipCount = 0; xy < xymax; ++xy, ++skipCount) {
                    if (data[xy] != prev[xy]) {
                        break;
                    }
                }
                xy -= skipCount;

                // determine repeat count 
                short v = data[xy];
                for (repeatCount = 0; xy < xymax && repeatCount < 127; ++xy, ++repeatCount) {
                    if (data[xy] != v) {
                        break;
                    }
                }
                xy -= repeatCount;

                if (skipCount < 2 && xy + skipCount < xymax && repeatCount < 3) {
                    literalCount++;
                    if (literalCount > 126) {
                        out.write(literalCount); // Literal OP-code
                        out.writeShorts(data, xy - literalCount + 1, literalCount);
                        literalCount = 0;
                    }
                } else {
                    if (literalCount > 0) {
                        out.write(literalCount); // Literal OP-code
                        out.writeShorts(data, xy - literalCount, literalCount);
                        literalCount = 0;
                    }
                    if (xy + skipCount == xymax) {
                        // => we can skip until the end of the line without
                        //    having to write an op-code
                        xy += skipCount - 1;
                    } else if (skipCount >= repeatCount) {
                        while (skipCount > 254) {
                            out.write(0); // Skip OP-code
                            out.write(255);
                            xy += 254;
                            skipCount -= 254;
                        }
                        out.write(0); // Skip OP-code
                        out.write(skipCount + 1);
                        xy += skipCount - 1;
                    } else {
                        out.write(-repeatCount); // Repeat OP-code
                        out.writeShort(v);
                        xy += repeatCount - 1;
                    }
                }
            }

            // flush literal run
            if (literalCount > 0) {
                out.write(literalCount);
                out.writeShorts(data, xy - literalCount, literalCount);
                literalCount = 0;
            }

            out.write(-1);// End of line OP-code
        }


        // Complete the header
        long pos = out.getStreamPosition();
        out.seek(headerPos);
        out.writeInt((int) (pos - headerPos));
        out.seek(pos);
    }
    /** Encodes a 24-bit key frame.
     *
     * @param out The output stream. Must be set to Big-Endian.
     * @param data The image data.
     * @param offset The offset to the first pixel in the data array.
     * @param length The width of the image in data elements.
     * @param step The number to add to offset to get to the next scanline.
     */
    public void writeKey24(ImageOutputStream out, int[] data, int offset, int length, int step)
            throws IOException {
        long headerPos = out.getStreamPosition();

        // Reserve space for the header:
        out.writeInt(0);
        out.writeShort(0x0000);

        // Encode each scanline
        int ymax = data.length;
        for (int y = offset; y < ymax; y += step) {
            int xy = y;
            int xymax = y + length;

            out.write(1); // this is a key-frame, there is nothing to skip at the start of line

            int literalCount = 0;
            int repeatCount = 0;
            for (; xy < xymax; ++xy) {
                // determine repeat count
                int v = data[xy];
                for (repeatCount = 0; xy < xymax && repeatCount < 127; ++xy, ++repeatCount) {
                    if (data[xy] != v) {
                        break;
                    }
                }
                xy -= repeatCount;

                if (repeatCount < 3) {
                    literalCount++;
                    if (literalCount > 126) {
                        out.write(literalCount); // Literal OP-code
                        writeInts24(out,data, xy - literalCount + 1, literalCount);
                        literalCount = 0;
                    }
                } else {
                    if (literalCount > 0) {
                        out.write(literalCount); // Literal OP-code
                        writeInts24(out,data, xy - literalCount, literalCount);
                        literalCount = 0;
                    }
                    out.write(-repeatCount); // Repeat OP-code
                    writeInt24(out,v);
                    xy += repeatCount - 1;
                }
            }

            // flush literal run
            if (literalCount > 0) {
                out.write(literalCount);
                writeInts24(out,data, xy - literalCount, literalCount);
                literalCount = 0;
            }

            out.write(-1);// End of line OP-code
        }


        // Complete the header
        long pos = out.getStreamPosition();
        out.seek(headerPos);
        out.writeInt((int) (pos - headerPos));
        out.seek(pos);
    }

    /** Encodes a 24-bit delta frame.
     *
     * @param out The output stream. Must be set to Big-Endian.
     * @param data The image data.
     * @param prev The image data of the previous frame.
     * @param offset The offset to the first pixel in the data array.
     * @param length The width of the image in data elements.
     * @param step The number to add to offset to get to the next scanline.
     */
    public void writeDelta24(ImageOutputStream out, int[] data, int[] prev, int offset, int length, int step)
            throws IOException {

        int height = (data.length - offset) / step;

        // Determine whether we can skip lines at the beginning
        int ymin;
        int ymax = offset + height * step;
        scanline:
        for (ymin = offset; ymin < ymax; ymin += step) {
            int xy = ymin;
            int xymax = ymin + length;
            for (; xy < xymax; ++xy) {
                if (data[xy] != prev[xy]) {
                    break scanline;
                }
            }
        }


        if (ymin == ymax) {
            // => Frame is identical to previous one
            out.writeInt(4);
            return;
        }

        // Determine whether we can skip lines at the end
        scanline:
        for (; ymax > ymin; ymax -= step) {
            int xy = ymax - step;
            int xymax = ymax - step + length;
            for (; xy < xymax; ++xy) {
                if (data[xy] != prev[xy]) {
                    break scanline;
                }
            }
        }
        //System.out.println("AppleRLEEncoder ymin:" + ymin / step + " ymax" + ymax / step);

        // Reserve space for the header
        long headerPos = out.getStreamPosition();
        out.writeInt(0);

        if (ymin == offset && ymax == offset + height * step) {
            // => we can't skip any lines:
            out.writeShort(0x0000);
        } else {
            // => we can skip lines:
            out.writeShort(0x0008);
            out.writeShort(ymin / step);
            out.writeShort(0);
            out.writeShort((ymax - ymin + 1) / step);
            out.writeShort(0);
        }


        // Encode each scanline
        for (int y = ymin; y < ymax; y += step) {
            int xy = y;
            int xymax = y + length;

            // determine skip count
            int skipCount = 0;
            for (; xy < xymax; ++xy, ++skipCount) {
                if (data[xy] != prev[xy]) {
                    break;
                }
            }
            if (skipCount == length) {
                // => the entire line can be skipped
                out.write(1); // don't skip any pixels
                out.write(-1); // end of line
                continue;
            }
            out.write(Math.min(255, skipCount + 1));
            if (skipCount > 254) {
                skipCount -= 254;
                while (skipCount > 254) {
                    out.write(0); // Skip OP-code
                    out.write(255);
                    skipCount -= 254;
                }
                out.write(0); // Skip OP-code
                out.write(skipCount + 1);
            }

            int literalCount = 0;
            int repeatCount = 0;
            for (; xy < xymax; ++xy) {
                // determine skip count
                for (skipCount = 0; xy < xymax; ++xy, ++skipCount) {
                    if (data[xy] != prev[xy]) {
                        break;
                    }
                }
                xy -= skipCount;

                // determine repeat count
                int v = data[xy];
                for (repeatCount = 0; xy < xymax && repeatCount < 127; ++xy, ++repeatCount) {
                    if (data[xy] != v) {
                        break;
                    }
                }
                xy -= repeatCount;

                if (skipCount < 2 && xy + skipCount < xymax && repeatCount < 3) {
                    literalCount++;
                    if (literalCount > 126) {
                        ImageOutputStreamImpl impl;
                        out.write(literalCount); // Literal OP-code
                        writeInts24(out,data, xy - literalCount + 1, literalCount);
                        literalCount = 0;
                    }
                } else {
                    if (literalCount > 0) {
                        out.write(literalCount); // Literal OP-code
                        writeInts24(out,data, xy - literalCount, literalCount);
                        literalCount = 0;
                    }
                    if (xy + skipCount == xymax) {
                        // => we can skip until the end of the line without
                        //    having to write an op-code
                        xy += skipCount - 1;
                    } else if (skipCount >= repeatCount) {
                        while (skipCount > 254) {
                            out.write(0); // Skip OP-code
                            out.write(255);
                            xy += 254;
                            skipCount -= 254;
                        }
                        out.write(0); // Skip OP-code
                        out.write(skipCount + 1);
                        xy += skipCount - 1;
                    } else {
                        out.write(-repeatCount); // Repeat OP-code
                        writeInt24(out,v);
                        xy += repeatCount - 1;
                    }
                }
            }

            // flush literal run
            if (literalCount > 0) {
                out.write(literalCount);
                writeInts24(out,data, xy - literalCount, literalCount);
                literalCount = 0;
            }

            out.write(-1);// End of line OP-code
        }


        // Complete the header
        long pos = out.getStreamPosition();
        out.seek(headerPos);
        out.writeInt((int) (pos - headerPos));
        out.seek(pos);
    }
    /** Encodes a 32-bit key frame.
     *
     * @param out The output stream. Must be set to Big-Endian.
     * @param data The image data.
     * @param offset The offset to the first pixel in the data array.
     * @param length The width of the image in data elements.
     * @param step The number to add to offset to get to the next scanline.
     */
    public void writeKey32(ImageOutputStream out, int[] data, int offset, int length, int step)
            throws IOException {
        long headerPos = out.getStreamPosition();

        // Reserve space for the header:
        out.writeInt(0);
        out.writeShort(0x0000);

        // Encode each scanline
        int ymax = data.length;
        for (int y = offset; y < ymax; y += step) {
            int xy = y;
            int xymax = y + length;

            out.write(1); // this is a key-frame, there is nothing to skip at the start of line

            int literalCount = 0;
            int repeatCount = 0;
            for (; xy < xymax; ++xy) {
                // determine repeat count
                int v = data[xy];
                for (repeatCount = 0; xy < xymax && repeatCount < 127; ++xy, ++repeatCount) {
                    if (data[xy] != v) {
                        break;
                    }
                }
                xy -= repeatCount;

                if (repeatCount < 3) {
                    literalCount++;
                    if (literalCount > 126) {
                        out.write(literalCount); // Literal OP-code
                        out.writeInts(data, xy - literalCount + 1, literalCount);
                        literalCount = 0;
                    }
                } else {
                    if (literalCount > 0) {
                        out.write(literalCount); // Literal OP-code
                        out.writeInts(data, xy - literalCount, literalCount);
                        literalCount = 0;
                    }
                    out.write(-repeatCount); // Repeat OP-code
                    out.writeInt(v);
                    xy += repeatCount - 1;
                }
            }

            // flush literal run
            if (literalCount > 0) {
                out.write(literalCount);
                out.writeInts(data, xy - literalCount, literalCount);
                literalCount = 0;
            }

            out.write(-1);// End of line OP-code
        }


        // Complete the header
        long pos = out.getStreamPosition();
        out.seek(headerPos);
        out.writeInt((int) (pos - headerPos));
        out.seek(pos);
    }

    /** Encodes a 32-bit delta frame.
     *
     * @param out The output stream. Must be set to Big-Endian.
     * @param data The image data.
     * @param prev The image data of the previous frame.
     * @param offset The offset to the first pixel in the data array.
     * @param length The width of the image in data elements.
     * @param step The number to add to offset to get to the next scanline.
     */
    public void writeDelta32(ImageOutputStream out, int[] data, int[] prev, int offset, int length, int step)
            throws IOException {

        int height = (data.length - offset) / step;

        // Determine whether we can skip lines at the beginning
        int ymin;
        int ymax = offset + height * step;
        scanline:
        for (ymin = offset; ymin < ymax; ymin += step) {
            int xy = ymin;
            int xymax = ymin + length;
            for (; xy < xymax; ++xy) {
                if (data[xy] != prev[xy]) {
                    break scanline;
                }
            }
        }


        if (ymin == ymax) {
            // => Frame is identical to previous one
            out.writeInt(4);
            return;
        }

        // Determine whether we can skip lines at the end
        scanline:
        for (; ymax > ymin; ymax -= step) {
            int xy = ymax - step;
            int xymax = ymax - step + length;
            for (; xy < xymax; ++xy) {
                if (data[xy] != prev[xy]) {
                    break scanline;
                }
            }
        }
        //System.out.println("AppleRLEEncoder ymin:" + ymin / step + " ymax" + ymax / step);

        // Reserve space for the header
        long headerPos = out.getStreamPosition();
        out.writeInt(0);

        if (ymin == offset && ymax == offset + height * step) {
            // => we can't skip any lines:
            out.writeShort(0x0000);
        } else {
            // => we can skip lines:
            out.writeShort(0x0008);
            out.writeShort(ymin / step);
            out.writeShort(0);
            out.writeShort((ymax - ymin + 1) / step);
            out.writeShort(0);
        }


        // Encode each scanline
        for (int y = ymin; y < ymax; y += step) {
            int xy = y;
            int xymax = y + length;

            // determine skip count
            int skipCount = 0;
            for (; xy < xymax; ++xy, ++skipCount) {
                if (data[xy] != prev[xy]) {
                    break;
                }
            }
            if (skipCount == length) {
                // => the entire line can be skipped
                out.write(1); // don't skip any pixels
                out.write(-1); // end of line
                continue;
            }
            out.write(Math.min(255, skipCount + 1));
            if (skipCount > 254) {
                skipCount -= 254;
                while (skipCount > 254) {
                    out.write(0); // Skip OP-code
                    out.write(255);
                    skipCount -= 254;
                }
                out.write(0); // Skip OP-code
                out.write(skipCount + 1);
            }

            int literalCount = 0;
            int repeatCount = 0;
            for (; xy < xymax; ++xy) {
                // determine skip count
                for (skipCount = 0; xy < xymax; ++xy, ++skipCount) {
                    if (data[xy] != prev[xy]) {
                        break;
                    }
                }
                xy -= skipCount;

                // determine repeat count
                int v = data[xy];
                for (repeatCount = 0; xy < xymax && repeatCount < 127; ++xy, ++repeatCount) {
                    if (data[xy] != v) {
                        break;
                    }
                }
                xy -= repeatCount;

                if (skipCount < 2 && xy + skipCount < xymax && repeatCount < 3) {
                    literalCount++;
                    if (literalCount > 126) {
                        out.write(literalCount); // Literal OP-code
                        out.writeInts(data, xy - literalCount + 1, literalCount);
                        literalCount = 0;
                    }
                } else {
                    if (literalCount > 0) {
                        out.write(literalCount); // Literal OP-code
                        out.writeInts(data, xy - literalCount, literalCount);
                        literalCount = 0;
                    }
                    if (xy + skipCount == xymax) {
                        // => we can skip until the end of the line without
                        //    having to write an op-code
                        xy += skipCount - 1;
                    } else if (skipCount >= repeatCount) {
                        while (skipCount > 254) {
                            out.write(0); // Skip OP-code
                            out.write(255);
                            xy += 254;
                            skipCount -= 254;
                        }
                        out.write(0); // Skip OP-code
                        out.write(skipCount + 1);
                        xy += skipCount - 1;
                    } else {
                        out.write(-repeatCount); // Repeat OP-code
                        out.writeInt(v);
                        xy += repeatCount - 1;
                    }
                }
            }

            // flush literal run
            if (literalCount > 0) {
                out.write(literalCount);
                out.writeInts(data, xy - literalCount, literalCount);
                literalCount = 0;
            }

            out.write(-1);// End of line OP-code
        }


        // Complete the header
        long pos = out.getStreamPosition();
        out.seek(headerPos);
        out.writeInt((int) (pos - headerPos));
        out.seek(pos);
    }

    private byte[] byteBuf = new byte[3];
    private void writeInt24(ImageOutputStream out, int v) throws IOException {
        if (out.getByteOrder() == ByteOrder.BIG_ENDIAN) {
            byteBuf[0] = (byte)(v >>> 16);
            byteBuf[1] = (byte)(v >>>  8);
            byteBuf[2] = (byte)(v >>>  0);
        } else {
            byteBuf[0] = (byte)(v >>>  0);
            byteBuf[1] = (byte)(v >>>  8);
            byteBuf[2] = (byte)(v >>> 16);
        }
        out.write(byteBuf, 0, 3);
    }

    private void writeInts24(ImageOutputStream out, int[] i, int off, int len) throws IOException {
        // Fix 4430357 - if off + len < 0, overflow occurred
        if (off < 0 || len < 0 || off + len > i.length || off + len < 0) {
            throw new IndexOutOfBoundsException
                ("off < 0 || len < 0 || off + len > i.length!");
        }

        byte[] b = new byte[len*3];
        int boff = 0;
        if (out.getByteOrder() == ByteOrder.BIG_ENDIAN) {
            for (int j = 0; j < len; j++) {
                int v = i[off + j];
                //b[boff++] = (byte)(v >>> 24);
                b[boff++] = (byte)(v >>> 16);
                b[boff++] = (byte)(v >>> 8);
                b[boff++] = (byte)(v >>> 0);
            }
        } else {
            for (int j = 0; j < len; j++) {
                int v = i[off + j];
                //b[boff++] = (byte)(v >>> 0);
                b[boff++] = (byte)(v >>> 8);
                b[boff++] = (byte)(v >>> 16);
                b[boff++] = (byte)(v >>> 24);
            }
        }

        out.write(b, 0, len*3);
    }

    /*
    public static void main(String[] args) {
        short[] data = {//
            8, 1, 1, 1, 1, 2, 8,//
            8, 0, 2, 0, 0, 0, 8,//
            8, 2, 3, 4, 4, 3, 8,//
            8, 2, 2, 3, 4, 4, 8,//
            8, 1, 4, 4, 4, 5, 8};
        short[] prev = {//
            8, 1, 1, 1, 1, 1, 8, //
            8, 5, 5, 5, 5, 0, 8,//
            8, 3, 3, 3, 3, 3, 8,//
            8, 2, 2, 0, 0, 0, 8,//
            8, 2, 0, 0, 0, 5, 8};
        ByteArrayOutputStream buf = new ByteArrayOutputStream();
        ImageOutputStream out = new MemoryCacheImageOutputStream(buf);
        AppleRLEEncoder enc = new AppleRLEEncoder();
        try {
            enc.writeDelta16(out, data, prev, 1, 5, 7);
//            enc.writeKey16(out, data, 1, 5, 7);
            out.close();
            byte[] result = buf.toByteArray();
            System.out.println("size:" + result.length);
            System.out.println(Arrays.toString(result));

            System.out.print("0x [");
            for (int i = 0; i < result.length; i++) {
                if (i != 0) {
                    System.out.print(',');
                }
                String hex = "00" + Integer.toHexString(result[i]);
                System.out.print(hex.substring(hex.length() - 2));
            }
            System.out.println(']');
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }*/
}
