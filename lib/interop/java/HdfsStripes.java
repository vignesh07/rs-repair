// Writes one random stripe per code, encoded by Hadoop's own RS coder classes (RSUtil + GF256, unmodified),
// called the same way org.apache.hadoop.io.erasurecode.rawcoder.RSRawEncoder calls them.
// usage: java HdfsStripes <outdir> <shardlen> n k [n k ...]
import org.apache.hadoop.io.erasurecode.rawcoder.util.RSUtil;
import java.io.FileOutputStream;
import java.util.Random;

public class HdfsStripes {
    public static void main(String[] a) throws Exception {
        int len = Integer.parseInt(a[1]);
        Random rng = new Random(20261126L);
        for (int i = 2; i + 1 < a.length; i += 2) {
            int n = Integer.parseInt(a[i]), k = Integer.parseInt(a[i + 1]);
            byte[] encodeMatrix = new byte[n * k];
            RSUtil.genCauchyMatrix(encodeMatrix, n, k);
            byte[] gfTables = new byte[n * k * 32];
            RSUtil.initTables(k, n - k, encodeMatrix, k * k, gfTables);
            byte[][] in = new byte[k][len], out = new byte[n - k][len];
            for (byte[] s : in) rng.nextBytes(s);
            RSUtil.encodeData(gfTables, len, in, new int[k], out, new int[n - k]);
            try (FileOutputStream f = new FileOutputStream(a[0] + "/stripe_" + n + "_" + k + ".bin")) {
                for (byte[] s : in) f.write(s);
                for (byte[] s : out) f.write(s);
            }
        }
    }
}
