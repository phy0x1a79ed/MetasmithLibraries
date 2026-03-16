process {
	withName: 'p01__downloadBuscoLineage' {
		cpus = 1
		memory = { (2**(task.attempt-1)) * ('4.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('2hours' as Duration) }
	}
	withName: 'p02__merge_bams' {
		cpus = 8
		memory = { (2**(task.attempt-1)) * ('16.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('2hours' as Duration) }
	}
	withName: 'p03__braker3' {
		cpus = 16
		memory = { (2**(task.attempt-1)) * ('64.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('2days' as Duration) }
	}
	withName: 'p04__stringtie_merge' {
		cpus = 2
		memory = { (2**(task.attempt-1)) * ('4.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('1hours' as Duration) }
	}
	withName: 'p05__gffread_proteins' {
		cpus = 2
		memory = { (2**(task.attempt-1)) * ('4.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('1hours' as Duration) }
	}
	withName: 'p06__stringtie_quant' {
		cpus = 4
		memory = { (2**(task.attempt-1)) * ('4.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('1hours' as Duration) }
	}
	withName: 'p07__busco' {
		cpus = 8
		memory = { (2**(task.attempt-1)) * ('16.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('4hours' as Duration) }
	}
	withName: 'p08__eggnog_mapper' {
		cpus = 16
		memory = { (2**(task.attempt-1)) * ('32.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('1day' as Duration) }
	}
	withName: 'p09__pydeseq2' {
		cpus = 2
		memory = { (2**(task.attempt-1)) * ('8.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('1hours' as Duration) }
	}
	withName: 'p10__stringtie_count_matrix' {
		cpus = 2
		memory = { (2**(task.attempt-1)) * ('4.00 GB' as MemoryUnit) }
		time = { (2**(task.attempt-1)) * ('1hours' as Duration) }
	}
}
