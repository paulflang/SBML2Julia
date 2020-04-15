	def get_gif_from_inchi(self, alphabet):
		with open('dna_stuctures.html', 'w') as file:
			file.writelines(['<!DOCTYPE html>\n', '<html>\n', '<body>\n'])

			for item in alphabet:
				this_inchi = item['structure']
				line = '<img src="https://cactus.nci.nih.gov/chemical/structure/{}/image">\n'.format(this_inchi)
				file.write(line)
				file.write(item["id"])

			file.write('</body>\n')
			file.write('</html>')