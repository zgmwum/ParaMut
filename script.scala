// $ADAM_HOME/bin/adam-shell --master $SPARK_CLUSTER_HOST --driver-memory 10G --packages datastax:spark-cassandra-connector:2.4.0-s_2.11 --conf spark.cassandra.input.reads_per_sec=1000 --conf spark.cassandra.concurrent.reads=10

val file = "Blast_finished_domains.csv.bz2"
val blast = spark.read.format("csv").option("header", "true").option("inferSchema", true).option("delimiter", ",").load(file).select("gene_name","tx_id","AminoAcid_Pos","AminoAcid","para_id","para_gene","Para_AminoAcid_Pos","para_start_whole_domain","para_end_whole_domain","para_whole_sequence","start_whole_domain","end_whole_domain","protein_sequence").withColumnRenamed("AminoAcid_Pos","main_aapos").withColumnRenamed("Para_AminoAcid_Pos","orto_aapos").withColumnRenamed("tx_id","main_tx_id").withColumnRenamed("para_id","orto_tx_id")

val exome = spark.read.parquet("exome-vepped").withColumn("alleles", split($"allele_string","/")).withColumn("ref",$"alleles".getItem(0)).withColumn("alt",$"alleles".getItem(1)).drop("alleles")
val exomevep = exome.withColumn("x",explode($"transcript_consequences")).withColumn("chr", regexp_replace($"seq_region_name", "^chr", "")).select($"chr",$"start" as "pos",$"x.transcript_id" as "tid",$"x.gene_symbol" as "gene",$"x.hgvsp" as "hgvsp",$"x.protein_start" as "aapos", $"ref", $"alt").filter($"aapos".isNotNull).dropDuplicates()

val blastAndMainPos = blast.join(exomevep
    .withColumnRenamed("chr","srcchr")
    .withColumnRenamed("pos","srcpos")
    .withColumnRenamed("ref","srcref")
    .withColumnRenamed("alt","srcalt")
    .withColumnRenamed("gene","srcgene")
    .withColumnRenamed("hgvsp","srchgvsp")
    .withColumnRenamed("tid","main_tx_id")
    .withColumnRenamed("aapos","main_aapos"), Seq("main_tx_id","main_aapos"), "inner")
val blastAndSrcAndOrtoPos = blastAndMainPos.join(exomevep
    .withColumnRenamed("chr","ortchr")
    .withColumnRenamed("pos","ortpos")
    .withColumnRenamed("ref","ortref")
    .withColumnRenamed("alt","ortalt")
    .withColumnRenamed("gene","ortgene")
    .withColumnRenamed("hgvsp","orthgvsp")
    .withColumnRenamed("tid","orto_tx_id")
    .withColumnRenamed("aapos","orto_aapos"), Seq("orto_tx_id","orto_aapos"), "inner")

blastAndSrcAndOrtoPos.write.parquet("/ortomut-with-positions-refalt-final/prq")

import org.bdgenomics.adam.rdd.ADAMContext; sc.hadoopConfiguration.setBoolean("org.bdgenomics.adam.converters.VariantContextConverter.NEST_ANN_IN_GENOTYPES", true);  val ac=new ADAMContext(sc)

val blastAndSrcAndOrtoPos = spark.read.parquet("/mnt/zgmvol/_forge/PiotrS/vmurcia/homologi-do-annotacio/ortomut-with-positions-refalt-final/prq").filter((abs($"srcpos" - $"ortpos") > 30) || ($"srcchr"!==$"ortchr"))


blastAndSrcAndOrtoPos.groupBy("main_tx_id","main_aapos","orto_tx_id","orto_aapos").agg(min($"gene_name") , min($"AminoAcid") , min($"para_gene") , min($"para_start_whole_domain") , min($"para_end_whole_domain") , min($"para_whole_sequence") , min($"start_whole_domain") , min($"end_whole_domain") , min($"protein_sequence") , min($"srcchr") , min($"srcpos") , min($"srcgene")  ,  min($"ortchr") , min($"ortpos") , min($"ortgene")  ).repartition(1).write.option("header", "true").option("delimiter", "\t").csv("/mnt/zgmvol/_forge/vmurcia/workspace/2019_02_18-wyniki/2019-10-17/paramut.csv")

val cln = ac.loadVcf("/clinvar.vcf.gz")
val clnReady = cln.toVariants.toDF.
    filter($"annotation.attributes.CLNSIG" === "Pathogenic" ||$"annotation.attributes.CLNSIG" === "Likely_pathogenic" || $"annotation.attributes.CLNSIG" === "Uncertain_significance").
    withColumn("disease_cln",$"annotation.attributes.CLNDN").
    withColumn("clnId",$"annotation.attributes.RS").
    select($"annotation.attributes.CLNSIG" as "CLINSIG", $"disease_cln",$"clnId",$"referenceName" as "chr",($"start" + 1) as "pos",$"referenceAllele" as "ref",$"alternateAllele" as "alt")

val clinvarOrtomuted = clnReady.join(blastAndSrcAndOrtoPos
    .withColumnRenamed("srcref","ref")
    .withColumnRenamed("srcalt","alt")
    .withColumnRenamed("srcchr","chr")
    .withColumnRenamed("srcpos","pos"), Seq("chr","pos","ref","alt"), "inner")

val clinvarOrtomutedReady = clinvarOrtomuted.withColumn("source", concat_ws(":", $"chr",$"pos",$"ref",$"alt") ).drop("chr","pos","ref","alt").withColumnRenamed("ortchr","chr").withColumnRenamed("ortpos","pos").select("chr","pos","disease_cln","clnId","source").dropDuplicates().repartition(1).write.option("header", "true").option("delimiter", "\t").csv("/mnt/zgmvol/_forge/vmurcia/workspace/2019_02_18-wyniki/2019-10-15/clinvarbisCut.csv")

clinvarOrtomuted.drop("orthgvsp","srchgvsp").dropDuplicates().repartition(1).write.option("header", "true").option("delimiter", "\t").csv("/mnt/zgmvol/_forge/vmurcia/workspace/2019_02_18-wyniki/2019-10-15/clinvarbis.csv")

clinvarOrtomuted.drop("orthgvsp","srchgvsp").dropDuplicates().join(cln.toVariants.toDF.select($"referenceName" as "ortchr",($"start" + 1) as "ortpos"),Seq("ortchr","ortpos"),"leftanti").repartition(1).write.option("header", "true").option("delimiter", "\t").csv("clinvarbisWoClinvar.csv")

