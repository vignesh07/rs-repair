// Stand-in for Hadoop's annotation class, so that the two original coder files compile on their own.
package org.apache.hadoop.classification;

public class InterfaceAudience {
    public @interface Private {}
    public @interface Public {}
    public @interface LimitedPrivate { String[] value(); }
}
