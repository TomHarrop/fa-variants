import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._

class ScatterBaseRecalibrator extends QScript {
  @Input(doc="The reference file for the bam files.", shortName="R")
  var referenceFile: File = _

  @Input(doc="Bam files to genotype.", shortName="I")
  var bamFiles: Seq[File] = _

  @Input(doc="A database of known polymorphic sites", shortName="knownSites")
  var knownSites: Seq[File] = _

  @Output(doc="output", shortName="out")
  var outFile: File = _

  @Argument(doc="One or more genomic intervals over which to operate", shortName="L", required=false)
  var intervals: File = _

  @Argument(doc="Input covariates table file for on-the-fly base quality score recalibration", shortName="BQSR", required=false)
  var bqsr: File = _

  def script() {
    val br = new BaseRecalibrator
    // Run options
    br.read_buffer_size = 1000000
    br.nct = 2
    br.memoryLimit = 6
    br.javaGCThreads = 1
    br.scatterCount = 4
    // HaplotypeCaller options
    br.interval_padding = if (intervals == null) 0 else 100
    // Files
    br.R = referenceFile
    br.I = bamFiles
    br.knownSites = knownSites
    br.out = outFile
    br.BQSR = bqsr
    br.intervals = if (intervals == null) Nil else List(intervals)
    // Run HaplotypeCaller
    add(br)
  }
}
