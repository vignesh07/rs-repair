// Writes one random stripe per code, encoded by Backblaze's JavaReedSolomon exactly as shipped.
// usage: java BBStripes <outdir> <shardlen> n k [n k ...]
import com.backblaze.erasure.ReedSolomon;
import java.io.FileOutputStream;
import java.util.Random;

public class BBStripes {
    public static void main(String[] a) throws Exception {
        int len = Integer.parseInt(a[1]);
        Random rng = new Random(20261126L);
        for (int i = 2; i + 1 < a.length; i += 2) {
            int n = Integer.parseInt(a[i]), k = Integer.parseInt(a[i + 1]);
            byte[][] shards = new byte[n][len];
            for (int s = 0; s < k; s++) rng.nextBytes(shards[s]);
            ReedSolomon rs = ReedSolomon.create(k, n - k);
            rs.encodeParity(shards, 0, len);
            if (!rs.isParityCorrect(shards, 0, len)) throw new RuntimeException("parity check failed");
            try (FileOutputStream f = new FileOutputStream(a[0] + "/stripe_" + n + "_" + k + ".bin")) {
                for (byte[] s : shards) f.write(s);
            }
        }
    }
}
