

import torch
import matplotlib.pyplot as plt #for plotting

from torchvision import datasets, transforms




def setup_device():
	torch.cuda.set_device(0)
	print('Used GPU Name:', torch.cuda.get_device_name(torch.cuda.current_device()))


def load_images(path):
	mTransform = transforms.Compose([transforms.Resize(512), transforms.ToTensor()])
	
	dataset = datasets.ImageFolder(root=path, transform=mTransform)
	
	dataloader = torch.utils.data.DataLoader(dataset, batch_size=16, shuffle=False)

	
	images_bad, labels = next(iter(dataloader))
	print(images_bad.size())
	images = images_bad[:,0,:,:]
	print(images.size())
	
	plt.figure(figsize=(3.200, 2.400), dpi=1000)
	plt.imshow(images[6])

	
	
setup_device()
	
load_images('Images/')