import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode.GVCF

class ScatterHaplotypeCaller extends QScript {
  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _

  @Input(doc="Bam files to genotype.", shortName="I")
  var bamFiles: Seq[File] = _

  @Output(doc="VCF output", shortName="out")
  var outFile: File = _

  @Argument(doc="One or more genomic intervals over which to operate", shortName="L", required=false)
  var intervals: File = _

  def script() {
    val hc = new HaplotypeCaller
    // Run options
    hc.read_buffer_size = 1000000
    hc.nct = 2
    hc.memoryLimit = 6
    hc.javaGCThreads = 1
    hc.scatterCount = 4
    // HaplotypeCaller options
    hc.interval_padding = if (intervals == null) 0 else 100
    hc.dontUseSoftClippedBases = true
    hc.emitRefConfidence = GVCF
    hc.stand_call_conf = 20.0
    hc.stand_emit_conf = 20.0
    // Files
    hc.R = referenceFile
    hc.I = bamFiles
    hc.out = outFile
    hc.intervals = if (intervals == null) Nil else List(intervals)
    // Run HaplotypeCaller
    add(hc)
  }
}
